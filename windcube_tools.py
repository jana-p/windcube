import datetime as dt
import time

import numpy as np
from scipy import optimize
import pandas as pd
import pdb
import multiprocessing as mp

# contains windcube constants and custom settings
import config_lidar as cl
import windcube_io as wio
import windcube_plotting as wp



# changes single LOS scan IDs of VAD composites to one ID
def change_scan_IDs(df):
    df['new_scan_ID'] = df['scan_ID'].copy()
    for cScan in cl.ScanID['COM']:
        df.loc[df['scan_ID'].isin( cl.CompDict[cScan] ), 'new_scan_ID'] = cScan

    df.drop(['scan_ID'], axis=1, inplace=True)
    df.rename(columns={'new_scan_ID':'scan_ID'}, inplace=True)

    return df


def run_fit(wrbin, az, ele, sProp, rbin):
    clms = ['wspeed', 'w', 'wdir', 'number_of_function_calls', 'rsquared']
    ws_out = pd.DataFrame( columns=clms, index=[rbin], dtype='float64' )
    # Target function
    fitfunc = lambda p, x: p[0]+p[1]*np.cos(x-p[2])
    # p[0] ... a (offset)
    # p[1] ... b (amplitude)
    # p[2] ... theta_max (phase shift)
    # x ... azimuth angle (theta)
    # y ... radial wind 
    # Distance to target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    # set azimuth to range from 0 to 360 instead of 0 to -0
    theta = np.radians( az[az < max(az)] )
    elevation = np.radians( ele[0] )
    # fit originally for radial wind positive towards lidar, radial wind 
    # however changed in get_data to positive away from lidar
    # radial wind changed back here to use fit as it is
    radial_wind = wrbin[az < max(az)] * (-1.0)
#   we.set_outliers_to_nan(radial_wind)
    # initial guess
    if radial_wind.any():
        # Initial guess for fit parameters
        guess_a = np.median( radial_wind )
        guess_b = 3*np.std( radial_wind )/(2**0.5)
        guess_phase = theta[ radial_wind.argmax(axis=0) ]
        p0 = [guess_a, guess_b, guess_phase]
        # least square fit
        try:
            p1, C, info, mes, success = optimize.leastsq(errfunc, p0[:], \
                    args=(theta, radial_wind), full_output=1)
            # calculate R^2
            ss_err=(info['fvec']**2).sum()
            ss_tot=((radial_wind-radial_wind.mean())**2).sum()
            rsquared=1-(ss_err/ss_tot)
            # number of function calls
            nfcalls = info['nfev']
        except ValueError:
            wio.printif( '..... ValueError in leastsq fit function' )
            rsquared = -999
            nfcalls = -999
        except TypeError:
            wio.printif( '..... TypeError in leastsq fit function' )
            rsquared = -999
            nfcalls = -999
        except RuntimeWarning:
            wio.printif( '..... nan in fit results' )
            rsquared = -999
            nfcalls = -999
    else:
        rsquared = -999
        nfcalls = -999
    ws_out['number_of_function_calls'][rbin] = nfcalls
    ws_out['rsquared'][rbin] = rsquared

    if rsquared>=0.1:
        # wind components
        # horizontal wind
        ws_out['wspeed'][rbin] = p1[1]/np.cos( elevation )
        # vertical wind
        ws_out['w'][rbin] = -p1[0]/np.sin( elevation )
        # wind direction
        ws_out['wdir'][rbin] = np.degrees(p1[2])
    else:
        ws_out['wspeed'][rbin] = np.nan
        ws_out['w'][rbin] = np.nan
        ws_out['wdir'][rbin] = np.nan

    return ws_out


# run loop over all VAD scans and fits sine function at all ranges
def wind_fit(AllW, sProp, sDate):
    if cl.SWITCH_POOL>0 and cl.SWITCH_POOL<5:
        # open number of pools specified in config file (SWITCH_POOL)
        wio.printif( '.... open small pool ' )
        pool = mp.Pool(processes=cl.SWITCH_POOL)
        POOL_LARGE = False
    elif cl.SWITCH_POOL>5:
        # open number of pools specified in config file (SWITCH_POOL)
        wio.printif( '.... open large pool ' )
        pool = mp.Pool(processes=cl.SWITCH_POOL)
        POOL_LARGE = True
    else:
        pool = 'dummy'
        POOL_LARGE = False

    # runs per range bin for large pool, or per VADscan for small pool
    scanIDs = AllW['scan_ID'].unique()
    if POOL_LARGE or pool=='dummy' or len(scanIDs)==1:
        combodf = pd.concat([ fit_parallel(AllW, sProp, sDate, VADscan, pool) \
                for VADscan in cl.ScanID['VAD'] ])
    else:
        poolres = [ pool.apply_async(fit_parallel,\
                args=(AllW, sProp, sDate, VADscan, 'dummy'))\
                for VADscan in cl.ScanID['VAD'] ]
        combodf = pd.concat([res.get() for res in poolres])

    if cl.SWITCH_POOL>0:
        wio.printif( '.... close pool ' )
        pool.close()
        pool.join()
    
    # combine VAD scans at 15 and 75 degrees elevation
    if cl.SWITCH_INPUT=='text':
        combodf = combodf[
                ((combodf['ele']==cl.LowerAngle) 
                    & (combodf.index.get_level_values('range')<=cl.CombiAlt)) \
                | ((combodf['ele']==cl.UpperAngle) 
                    & (combodf.index.get_level_values('range')>cl.CombiAlt))
                ]

        subs = ['w','wspeed','wdir']
        combodf.dropna( axis=0, how='all', subset=subs, inplace=True )
        combodf = combodf.reset_index()
        combodf['time'] = pd.to_datetime( combodf['time'], unit='s' )
        combodf = combodf.set_index(['time','range']).unstack(level='range').resample('15T').mean().stack(level='range')

        # plot timeseries of combined VAD
        if cl.SWITCH_PLOT:
            anglestr = str(cl.LowerAngle) + '-' + str(cl.UpperAngle)
            plotvars = [
                    ['wspeed','horizontal wind speed / m/s',anglestr],
                    ['w','vertical wind speed / m/s (positive = updraft)',anglestr],
                    ['wdir','wind direction / degrees (0, 360 = North)',anglestr],
                    ]
            for pv in plotvars:
                wp.plot_ts(combodf, sProp, sDate, pv)


def fit_parallel( AllW, sProp, sDate, VADscan, pool ):
    w = AllW[AllW.scan_ID==VADscan]
    if len(w.scan_ID)>0:
        wio.printif('.... fitting VAD ' + str(VADscan) )
        # separate different scans
        t=[x[0] for x in w.index]
        r=[x[1] for x in w.index]
        newScanIx=np.where(np.diff(t)>dt.timedelta(seconds=59))
        newScanIx = newScanIx[0]
        newScanPlus=np.concatenate([[0],newScanIx[0:-1]])
        # repeat fitting procedure for each VAD scan
        onerange = w.index.get_level_values('range')[
                w.index.get_level_values('time')==w.index.get_level_values('time')[0]]
        windindex = pd.MultiIndex.from_product(
                [w.index.get_level_values('time')[newScanPlus],onerange], \
                        names=['time', 'range'])
        clms = [ 'wspeed', 'w', 'wdir', 
                'confidence_index', 'rsquared', 'number_of_function_calls' ]
        win = pd.DataFrame(index=windindex, columns=clms, dtype='float64')
        if cl.SWITCH_POOL>0 and pool<>'dummy':
            # fit each VAD scan parallel in different pool
            poolres = [
                    pool.apply_async(loop_pool, \
                            args=(w, ixs, newScanIx, newScanPlus, win, sProp) )\
                            for ixs in range(0, len(newScanIx))
                            ]
            wind = pd.concat([res.get() for res in poolres])
        elif cl.SWITCH_POOL==0 or pool=='dummy':
            # run loop if number of pools is 0 (not parallel)
            windlist = [loop_pool(w, ixs, newScanIx, newScanPlus, win, sProp)\
                    for ixs in range(0, len(newScanIx)) ]
            wind = pd.concat(windlist)
        else:
            wio.printif( '... Please check SWITCH_POOL in config file!' )

        wind.set_index(windindex, inplace=True)
        # change negative wind direction (adapt speed as well)
        wind.loc[wind.wdir<0, 'wspeed'] = np.absolute( wind.loc[wind.wdir<0, 'wspeed'] )
        wind.loc[wind.wdir<0, 'wdir'] = np.absolute( wind.loc[ wind.wdir<0, 'wdir' ] )
        elestr = str( int( round( w['ele'][0] ) ) )
        # change range to altitude above ground level
        wind['alt'] = [x[1] * np.sin( np.radians( float(elestr) ) ) for x in wind.index]
        clms = {'range': 'oldrange', 'alt': 'range'}
        idx = ['time','range']
        wind = wind.reset_index().rename(columns=clms).set_index(idx)

        if cl.SWITCH_PLOT:
            plotvars = [
#                   ['rsquared', 'rsquared', elestr],
#                   ['confidence_index', 'confidence_index', elestr],
                    ['wspeed', 'horizontal wind speed / m/s', elestr, VADscan],
                    ['w', 'vertical wind speed / m/s (positive = updraft)', elestr, VADscan],
                    ['wdir', 'wind direction / degrees (0, 360 = North)', elestr, VADscan],
                    ]
            for pv in plotvars:
                wp.plot_ts(wind, sProp, sDate, pv)

        if cl.SWITCH_OUTNC:
            wind, pkllist = wio.get_pickle(sProp, ['VAD',elestr], wind)
            if not pkllist:
                wio.write_pickle(wind, [cl.ncInput + 'VAD' + '_' + elestr + '.pkl'], [elestr])
            else:
                wio.write_pickle(wind, [pkllist[0]], [elestr])
            wind = wind[~wind.index.duplicated(keep='last')]
            wio.export_to_netcdf(wind, sProp, sDate, '_VAD' + '_' + elestr)
                
        wind['ele'] = int( elestr )

        return wind


def loop_pool(w, ixs, newScanIx, newScanPlus, win, sProp):
    s0 = newScanPlus[ixs]
    s = newScanIx[ixs]
    ws = w[s0:s+1]
    meanconf = ws.confidence_index.mean(level=1)
    elevation = np.radians( ws['ele'][0] )
    # run fit for each height bin
    gr = ws.groupby(level='range')
    fitlist=[]
    gr.apply(lambda x: fitlist.append( run_fit(
        x[cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]], 
        x['azi'], x['ele'], sProp, 
        x.index.get_level_values('range')[0])) )

    wfit_out = pd.concat( fitlist )
    # filter out fits with rsquared smaller 0.5 using the confidence index
    wfit_out['confidence_index'] = meanconf
    wfit_out.loc[wfit_out['rsquared']<0.5, 'confidence_index'] = 0

    return wfit_out


# calculates time since "oldtime" and prints if output option is set to True
def timer(oldtime):
    if cl.SWITCH_TIMER:
        newtime = time.time()
        wio.printif(dt.timedelta(seconds=newtime - oldtime))
        return newtime

