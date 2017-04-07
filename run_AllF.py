import datetime as dt
import os
import pandas as pd
from glob import glob
import pdb
import multiprocessing as mp

# contains all windcube constants:
import config_lidar as cl

# contain windcube functions:
import windcube_io as wio
import windcube_tools as wt
import windcube_extras as we
import windcube_plotting as wp

if cl.SWITCH_TIMER:
    import time

#import warnings
#warnings.simplefilter("ignore", FutureWarning)
#warnings.simplefilter("error")



def main(sDate,p):
    STARTTIME = time.time()
    AllF = pd.DataFrame()

    # find files to read
    fend = cl.VarDict[p]['cols'][cl.VarDict[p]['N']]
    InNC = sorted(glob(cl.ncInput + fend + '*.nc'))
    InPKL = sorted(glob(cl.ncInput + fend + '.pkl'))
    InTXT = sorted(glob(cl.txtInput + cl.VarDict[p]['fend'] + '.' + cl.ending))

    # read input depending on settings in config file
    if InNC and (cl.SWITCH_INPUT=='append' or cl.SWITCH_INPUT=='netcdf'):
        # 1) read netcdf files in list
        wio.printif('... netcdf file found')
        AllF = pd.concat([AllF.append(wio.open_existing_nc(ncf)) \
                for ncf in InNC])
        # rename variables for consistent input
        if 'dv' in AllF:
            AllF.rename(columns={'dv':'radial_wind_speed'}, inplace=True)
        AllF.sort(inplace=True)
        # 2) appending txt to existing netcdf file
        if cl.SWITCH_INPUT=='append':
            wio.printif('... appending latest ascii file')
            AllF = AllF.append(wio.get_data(InTXT[-1],p))
            wio.printif('... sorting')
            AllF = AllF.sort_index()
            if cl.SWITCH_CLEANUP:
                os.remove(InTXT[-1])
        AllF.drop_duplicates(inplace=True)
    # 3) open pickle file
    elif (cl.SWITCH_INPUT=='text' or cl.SWITCH_INPUT=='pickle'):
        if InPKL and cl.SWITCH_INPUT=='pickle':
            wio.printif('... open pickle')
            AllF = pd.read_pickle(InPKL[0])
        # 4) read all available txt files
        # run parallel ... 
        if InTXT and cl.SWITCH_POOL>0:
            # open number of pools specified in config file (SWITCH_POOL)
            wio.printif( '... open pool ' )
            pool = mp.Pool(processes=cl.SWITCH_POOL)
            # read in files parallel
            poolres = [pool.apply_async(wio.get_data,args=(f,p)) for f in InTXT]
            pool.close()
            pool.join()
            wio.printif( '... close pool ' )
            # put results of parallel process in data frame
            AllFlist = [res.get() for res in poolres]
            AllF = AllF.append( pd.concat(AllFlist) )
            if cl.SWITCH_INPUT=='pickle':
                AllF.drop_duplicates(inplace=True)
                if cl.SWITCH_CLEANUP:
                    [os.remove(f) for f in InTXT]
        # ... or not parallel
        elif InTXT and cl.SWITCH_POOL==0:
            wio.printif('... loop over ascii files')
            # run loop if number of pools is 0 (not parallel)
            for f in InTXT:
                wio.printif('... reading file ' + f)
                AllF = AllF.append(wio.get_data(f,p))
                # remove hourly text files at the end of the day
                if ((len(InTXT)==24 and cl.SWITCH_OUTNC)\
                        or cl.SWITCH_INPUT=='pickle') and (cl.SWITCH_CLEANUP):
                    os.remove(f)
        elif not InTXT:
            wio.printif( '... no new text file' )
        else:
            wio.printif( '... Please check SWITCH_POOL in config file!' )
        
    INPUTTIME = wt.timer(STARTTIME)
#   pdb.set_trace()

    # change scan IDs of LOS to composite VAD
    AllF = wt.change_scan_IDs(AllF)

    # export content of data frame to netcdf
    if cl.SWITCH_OUTNC:
        wio.printif('... export nc')
        wio.export_to_netcdf(AllF,p,sDate,'')
        EXPORTTIME = wt.timer(STARTTIME)
    # export content of data frame to pickle file
    if cl.SWITCH_INPUT=='pickle':
        wio.printif('... export pickle')
        AllF.to_pickle(cl.ncInput + fend + '.pkl')
        # remove yesterday's pickle
        dtDate = dt.datetime.strptime( sDate, '%Y%m%d' ) - dt.timedelta(days=1)
        yesteryear = dt.datetime.strftime( dtDate, '%Y')
        yesterday = dt.datetime.strftime( dtDate, '%Y%m%d')
        pklfilename = os.path.split( os.path.split(cl.ncInput)[0] )[0] \
                + os.sep + yesteryear + os.sep + yesterday + '_' + fend + '.pkl'
        oldpkl = sorted(glob( pklfilename ))
        if oldpkl and cl.SWITCH_CLEANUP:
            wio.printif('... remove old pickle')
            [os.remove(f) for f in oldpkl]

    # plot time series (vertical line-of-sight only, 24h)
    if p<>'spectra':
        if ('LOS90' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                and cl.SWITCH_PLOT:
            if p=='wind':
                print('time series plot speed test')
                TSTIME1 = wt.timer(STARTTIME)
                wp.plot_ts(AllF,'cnr',sDate,['dummy'])
                TSTIME2 = wt.timer(TSTIME1)
                wp.plot_ts(AllF,'cnr',sDate,['dummy','fast'])
                TSTIME3 = wt.timer(TSTIME2)
            wp.plot_ts(AllF,p,sDate,['dummy'])
            TSTIME = wt.timer(STARTTIME)
    
        if p<>'dbs':
            # plot low level scans (polar)
            if ('LOW' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                    and cl.SWITCH_PLOT:
                [wp.plot_low_scan( AllF[AllF.scan_ID==LOWscan], p, sDate ) \
                        for LOWscan in cl.ScanID['LOW']]
                if p=='wind':
                    [wp.plot_low_scan( AllF[AllF.scan_ID==LOWscan], 'cnr', sDate ) \
                            for LOWscan in cl.ScanID['LOW']]
                LOWTIME = wt.timer(STARTTIME)

            # plot line-of-sight scans (scan duration)
            if ('LOS' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                    and cl.SWITCH_PLOT:
                [wp.plot_los( AllF[AllF.scan_ID==LOSscan], p, sDate ) \
                        for LOSscan in cl.ScanID['LOS']]
                if p=='wind':
                    [wp.plot_los( AllF[AllF.scan_ID==LOSscan], 'cnr', sDate ) \
                            for LOSscan in cl.ScanID['LOS']]
                LOSTIME = wt.timer(STARTTIME)

            # calculate horizontal wind speed and direction
            if 'VAD' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE:
                if p=='wind':
                    wio.printif('... fitting radial wind')
                    wt.wind_fit(AllF, p, sDate)

        if p=='dbs':
            # compare DBS wind components to VAD scan results 
            # (75 degrees elevation VAD scan only)
            we.compare_dbs(AllF, p, sDate)

    ENDTIME = wt.timer(STARTTIME)

    # create hdcp2-type output files if specified in config file 
    # (additional to standard netcdf output)
    if (cl.SWITCH_PLOT or cl.SWITCH_OUTNC) and cl.SWITCH_HDCP2:
        we.create_hdcp2_output(cl.sDate)



if __name__=="__main__":
    # creating output path if not existing
    if os.path.exists(cl.OutPath):
        if not os.path.exists(cl.ncOUT):
            os.makedirs(cl.ncOUT)
    else:
        print("Warning: Output path doesn't exist!")

    # running "main" function, starting processing
    for p in cl.proplist:
        if cl.SWITCH_PLOT or cl.SWITCH_OUTNC:
            main(cl.sDate,p)
            wio.printif('.. finished ' + p)
        else:
            print( '.. no output selected. Abort execution! \
                    Please adjust run options in config_lidar.py.' )

