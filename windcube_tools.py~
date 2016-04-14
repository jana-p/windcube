import datetime as dt
import time
import os

import numpy as np
from scipy import optimize
import pandas as pd
import pdb
import xray

import config_lidar as cl                               # contains all constants
import windcube_plotting as wp


# read text files from MySQL data base, returns pandas data frame
def get_data(file_path, sProp):
    # dsb text file have a different time stamp
    if sProp=='dbs':
        dparse = lambda d: dt.datetime.strptime(d, '%Y-%m-%d %H:%M:%S')
    else:
        dparse = lambda d: dt.datetime.strptime(d, '%Y-%m-%d %H:%M:%S.%f')
    
    outdf = pd.read_csv(file_path, sep='\t',             # read file from csv
            header=None, 
            skiprows=1,
            names=cl.VarDict[sProp]['cols'],
            index_col=['time', 'range'],                # index are time and range
            parse_dates=[0],                            # feed the first column to the parser
            date_parser=dparse, 
            squeeze=True                                # convert to `Series` object because we only have one column
            )
    # convert azimuth range from -180 - 180 degrees to 0 - 360 degrees
    outdf.loc[outdf.azi < 0, 'azi'] = outdf.loc[outdf.azi < 0, 'azi'] + 360
    # convert radial wind 'positive towards lidar' to radial wind 'positive away from lidar
    if 'radial_wind_speed' in outdf:
        outdf.radial_wind_speed = outdf.radial_wind_speed * (-1.0)

    return outdf


# opens existing netcdf and returns pandas data frame
def open_existing_nc(InFile):
    # open netcdf files
    xds = xray.open_dataset(InFile, decode_times=True ).transpose()
    # convert to pandas data frame
    dfnc = xds.to_dataframe()
    # swap order of indices (first time, then range)
    dfswap = dfnc.swaplevel('time','range')
    xds.close()

    return dfswap


# adds attribute to netcdf variable, if field is not empty
def add_attribute(varatts, v, where, what):
    if what<>'':
        varatts[where] = what

    return varatts


# adds all attributes to a dictionary, which is then added to the data frame
def all_att_to_df(df, df1D, vname, p, n):
    if vname in cl.AttDict:
        varatts = {}
        varatts = add_attribute( varatts, vname, 'units', cl.AttDict[vname][3] )
        varatts = add_attribute( varatts, vname, 'long_name', cl.AttDict[vname][2] )
        varatts = add_attribute( varatts, vname, 'standard_name', cl.AttDict[vname][1] )
        varatts = add_attribute( varatts, vname, 'comments', cl.AttDict[vname][5] )
        vnamenew = cl.AttDict[vname][0]
        df.rename( name_dict={ vname : vnamenew }, inplace=True )
        if cl.AttDict[vname][6]==1:
            df1D[vnamenew] = df[vnamenew][:,0]
            df1D[vnamenew].attrs = varatts
            df = df.drop( [vnamenew] )
        else:
            df[vnamenew].attrs = varatts
    else:
        if p<>'hdcp2':
            df[vname].attrs={
                    'units':cl.VarDict[p]['units'][n],
                    'long_name':cl.VarDict[p]['longs'][n], 
                    'standard_name':cl.VarDict[p]['names'][n]
                    }

    return df, df1D


# changes single LOS scan IDs of VAD composites to one ID
def change_scan_IDs(df):
    sID = df.scan_ID.unique()
    for cScan in cl.ScanID['COM']:
        df['scan_ID'][df['scan_ID'].isin( cl.CompDict[cScan] )] = cScan

    return df


# prepares pandas data frame for export to netcdf file
def export_to_netcdf(df,sProp,sDate,nameadd):
    # put spectra data in 3 dimensions (time, range, frequency)
    if sProp=='spectra':
        specS = df.spectra.str.split(',', expand=True).astype(float).stack()
        specdf = specS.to_frame()
        # add frequency index
        specdf.reset_index(inplace=True)
        specdf.rename( columns={ 'level_2' : 'frequency_bins' }, inplace=True )
        specdf.set_index(['time','range','frequency_bins'], inplace = True)
        fbix = specdf.unstack().columns.labels[1]
        vals = specdf.unstack().unstack().values
        vals = vals.reshape( np.shape(vals)[0], np.shape(vals)[1]/len(fbix), len(fbix) )
    printif('.... convert from df to xray ds')
    if nameadd == '':
        sID = df.scan_ID.unique()
        for s in sID:
            dfsID = df[df['scan_ID']==s]
            create_xray_dataset(dfsID, nameadd, s, sProp, sDate)
    elif 'VAD' in nameadd:
        create_xray_dataset(df, nameadd, 'VAD', sProp, sDate)


# exports xray data set to netcdf file, including global attributes, long names and units
def create_xray_dataset(df, nameadd, s, sProp, sDate):
    # change time index to seconds since 1970 for storing in netcdf
    printif('.... convert time to seconds since 1970')
    df.reset_index(inplace = True)
    tdt = df['time']
    t = df.time.astype(np.int64) / 10**9
    r = df['range']
    df['time'] = t
    df.set_index(['time','range'], inplace = True)
    xData = xray.Dataset.from_dataframe(df)
    xData1Ddict = {}
    # add variable attributes
    if 'VAD' in nameadd:
        pname='VAD'
    else:
        pname=sProp
    printif('.... add attributes to xray ds')
    if sProp=='hdcp2':
        for col in df:
            xData, xData1Ddict = all_att_to_df( xData, xData1Ddict, col, pname, [] )
    else:
        for c in range( 0, len(cl.VarDict[pname]['cols']) ):
            xData, xData1Ddict = all_att_to_df( xData, xData1Ddict, cl.VarDict[pname]['cols'][c], pname, c )
    # remove range coordinate from 1-dim data set and merge data sets
    if xData1Ddict=={}:
        xOut = xData
    else:
        xData1D = xray.Dataset(xData1Ddict).drop('range')
        xOut = xData.merge(xData1D)
        if sProp=='spectra':
            xData3D = xray.Dataset( { 'time' : ('time', xOut.time.values ),
                'range' : ('range', xOut['range'].values ),
                'frequency_bins' : ('frequency_bins', np.array(fbix) ),
                'spectra' : (['time','range','frequency_bins'], vals )
                } )
            xOut = xOut.drop( 'spectra' )
            xOut = xOut.merge(xData3D)
            xData3D.close()
            xOut[sProp].attrs={
                    'units':cl.VarDict[sProp]['units'][cl.VarDict[sProp]['N']],
                    'long_name':cl.VarDict[sProp]['longs'][cl.VarDict[sProp]['N']],
                    'standard_name':cl.VarDict[sProp]['names'][cl.VarDict[sProp]['N']]
                    }
        xData1D.close()
    xData.close()
    # add general variables (0-dim)
    for v in range( 0, len(cl.GenDict['cols']) ):
        xOut[cl.GenDict['cols'][v]] = cl.GenDict['val'][v]
        xOut[cl.GenDict['cols'][v]].attrs={
                'units':cl.GenDict['units'][v],
                'long_name':cl.GenDict['longs'][v],
                'standard_name':cl.GenDict['names'][v]
                }
    # add time stamp to dictionary of global attributes
    atts=cl.GloDict
    atts['Processing_date']=dt.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d, %H:%M:%S')
    # add global attributes to data set
    xOut.attrs=atts
    # set dummy date (short global attribute) to actual date
    xOut['year'] = np.int16( sDate[0:4] )
    xOut['month'] = np.int16( sDate[4:6] )
    xOut['day'] = np.int16( sDate[6:] )
#   pdb.set_trace()
    printif('.... write to file')
    # specify file name ending
    if sProp=='dbs':
        nameadd = 'DBS'
    elif sProp=='hdcp2':
        nameadd = 'HDCP2_' + nameadd
    else:
        if 'VAD' in nameadd:
            nameadd = nameadd[1:]
        else:
            nameadd = cl.VarDict[pname]['cols'][cl.VarDict[pname]['N']] + nameadd + '_scanID_' + str(s)
    # export file
    xOut.to_netcdf(path=cl.DataPath + sDate[0:4] + os.sep + sDate + '_' + nameadd + '.nc',
            mode='w', engine='netcdf4')#, format='NETCDF4')
    xOut.close()


# replaces outliers in data set with NaN (used before wind fit)
def set_outliers_to_nan(data_points):
    margin=40
    try:
        nd = np.abs(data_points - np.median(data_points))
        s = nd/np.median(nd)
        data_points[s>margin]=np.nan
    except IndexError:
        printif('.... wind fit: no outlier screening')

    return data_points


def run_fit(wrbin, az, ele, sProp, rbin):
    ws_out = pd.DataFrame( columns=['wspeed', 'w', 'wdir', 'number_of_function_calls', 'rsquared'], index=[rbin] )
    fitfunc = lambda p, x: p[0]+p[1]*np.cos(x-p[2]) # Target function
    # p[0] ... a (offset), p[1] ... b (amplitude), p[2] ... theta_max (phase shift)
    # x ... azimuth angle (theta)
    # y ... radial wind 
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    # set azimuth to range from 0 to 360 instead of 0 to -0
    theta = np.radians( az[az < max(az)] )
    elevation = np.radians( ele[0] )
    # fit originally for radial wind positive towards lidar, radial wind however changed in get_data to positive away from lidar
    # radial wind changed back here to use fit as it is
    radial_wind = wrbin[az < max(az)] * (-1.0)
    set_outliers_to_nan(radial_wind)
    # initial guess
    guess_a = np.median( radial_wind )
    guess_b = 3*np.std( radial_wind )/(2**0.5)
    guess_phase = theta[ radial_wind.argmax(axis=0) ]
    p0 = [guess_a, guess_b, guess_phase] # Initial guess for the parameters
    # least square fit
    try:
        p1, C, info, mes, success = optimize.leastsq(errfunc, p0[:], args=(theta, radial_wind), full_output=1)
        # calculate R^2
        ss_err=(info['fvec']**2).sum()
        ss_tot=((radial_wind-radial_wind.mean())**2).sum()
        rsquared=1-(ss_err/ss_tot)
        ws_out['number_of_function_calls'][rbin] = info['nfev'] # number of function calls
    except ValueError:
        rsquared=-999
        ws_out['number_of_function_calls'][rbin] = -999 # number of function calls
    ws_out['rsquared'][rbin] = rsquared # R^2

    if rsquared>=0.1:#0.2
        # wind components
        ws_out['wspeed'][rbin] = p1[1]/np.cos( elevation ) # horizontal wind
        ws_out['w'][rbin] = -p1[0]/np.sin( elevation ) # vertical wind
        ws_out['wdir'][rbin] = np.degrees(p1[2]) # wind direction
    else:
        # plot fit
#       pdb.set_trace()
#       plt.figure(figsize=(6, 5))
#       plt.scatter(theta,radial_wind)
#       plt.plot(theta, p1[0]+p1[1]*np.cos(theta-p1[2]))
#       plt.show()
#       plt.close()
        ws_out['wspeed'][rbin] = np.nan
        ws_out['w'][rbin] = np.nan
        ws_out['wdir'][rbin] = np.nan

    return ws_out


# run loop over all VAD scans and fits sine function at all ranges
def wind_fit(AllW, sProp, sDate):
    for VADscan in cl.ScanID['VAD']:
        printif('.... fitting VAD ' + str(VADscan) )
        # select only VAD scans
        w = AllW[AllW.scan_ID==VADscan]
        if len(w.scan_ID)>0:
            # separate different scans
            t=[x[0] for x in w.index]
            r=[x[1] for x in w.index]
            newScanIx=np.where(np.diff(t)>dt.timedelta(seconds=59))
            newScanPlus=np.concatenate([[0],newScanIx[0][0:-1]])
            # repeat fitting procedure for each VAD scan
            s0=0
            onerange = w.index.get_level_values('range')[w.index.get_level_values('time')==w.index.get_level_values('time')[0]]
            windindex = pd.MultiIndex.from_product([w.index.get_level_values('time')[newScanPlus],onerange], names=['time', 'range'])
            wind = pd.DataFrame(data=np.empty( [len(newScanPlus)*len(onerange), 6] ).fill(np.nan), 
                    index=windindex, 
                    columns=[ 'wspeed', 'w', 'wdir', 'confidence_index', 'rsquared', 'number_of_function_calls' ])
            for s in newScanIx[0]:
                ws = w[s0:s+1]
                meanconf = ws.confidence_index.mean(level=1)
                elevation = np.radians( ws['ele'][0] )
                # run fit for each height bin
                gr = ws.groupby(level='range')
                fitlist=[]
                gr.apply(lambda x: fitlist.append( run_fit(x[cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]], 
                        x['azi'], x['ele'], sProp, x.index.get_level_values('range')[0])) )
                wfit_out = pd.concat( fitlist )

                wind.loc[w.index.get_level_values('time')[s0], 'wspeed'] = wfit_out.wspeed.values
                wind.loc[w.index.get_level_values('time')[s0], 'w'] = wfit_out.w.values
                wind.loc[w.index.get_level_values('time')[s0], 'wdir'] = wfit_out.wdir.values
                wind.loc[w.index.get_level_values('time')[s0], 'rsquared'] = wfit_out.rsquared.values
                wind.loc[w.index.get_level_values('time')[s0], 'number_of_function_calls'] = wfit_out.number_of_function_calls.values
                wind.loc[w.index.get_level_values('time')[s0], 'confidence_index'] = meanconf.values
                s0 = s

            # change negative wind direction (adapt speed as well)
            wind.loc[ wind.wdir<0, 'wspeed' ] = np.absolute( wind.loc[ wind.wdir<0, 'wspeed' ] )
            wind.loc[ wind.wdir<0, 'wdir' ] = np.absolute( wind.loc[ wind.wdir<0, 'wdir' ] )
            if cl.SWITCH_PLOT:
                elestr = str( int( round( ws['ele'][0] ) ) )
                wp.plot_ts(wind, sProp, sDate, ['wspeed', 'horizontal wind speed / m/s', elestr])
                wp.plot_ts(wind, sProp, sDate, ['w', 'vertical wind speed / m/s (positive = updraft)', elestr])
                wp.plot_ts(wind, sProp, sDate, ['wdir', 'wind direction / degrees (0, 360 = North)', elestr])

            if cl.SWITCH_OUTNC:
                export_to_netcdf(wind, sProp, sDate, '_VAD' + '_' + elestr)


# prints message if output option is set in config file
def printif(message):
    if cl.SWITCH_OUTPUT:
        print(message)


# calculates time since "oldtime" and prints if output option is set to True
def timer(oldtime):
    if cl.SWITCH_TIMER:
        newtime = time.time()
        printif(dt.timedelta(seconds=newtime - oldtime))
        return newtime

