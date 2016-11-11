import datetime as dt
import os
import pandas as pd
from glob import glob
import pdb
import multiprocessing as mp

# contains all windcube constants:
import config_lidar as cl

# contain windcube functions:
import windcube_tools as wt
import windcube_extras as we
import windcube_plotting as wp

if cl.SWITCH_TIMER:
    import time


def main(sDate,p):
    STARTTIME = time.time()
    AllF = pd.DataFrame()
    # find files to read
    fend = cl.VarDict[p]['cols'][cl.VarDict[p]['N']]
    InNC = sorted(glob(cl.ncInput + fend + '.nc'))
    InTXT = sorted(glob(cl.txtInput + cl.VarDict[p]['fend'] + '.' + cl.ending))
    # read input depending on settings in config file
    # 1) appending txt to existing netcdf file
    if len(InNC)>=1 and cl.SWITCH_INNC=='append':
        wt.printif('... nc found')
        AllF = AllF.append(wt.open_existing_nc(InNC[0]))
        wt.printif('... appending latest ascii file')
        AllF = AllF.append(wt.get_data(InTXT[-1],p))
        wt.printif('... sorting')
        AllF = AllF.sort_index()
        os.remove(InTXT[-1])
    # 2) read (first) netcdf file in list
    elif len(InNC)>=1 and cl.SWITCH_INNC:
        AllF = AllF.append(wt.open_existing_nc(InNC[0]))
    # 3) read all available txt files
    else:
        wt.printif('... loop over ascii files')
        # run parallel ... 
        if cl.SWITCH_POOL>0:
            # open number of pools specified in config file (SWITCH_POOL)
            wt.printif( '... open pool ' )
            pool = mp.Pool(processes=cl.SWITCH_POOL)
            # read in files parallel
            poolres = [pool.apply_async(wt.get_data,args=(f,p)) for f in InTXT]
            pool.close()
            pool.join()
            wt.printif( '... close pool ' )
            counti = 0
            # put results of parallel process in data frame
            for res in poolres:
                if counti == 0:
                    AllF = res.get()
                else:
                    AllF = AllF.append( res.get() )
                counti = 1
        # ... or not
        elif cl.SWITCH_POOL==0:
            # run loop if number of pools is 0 (not parallel)
            for f in InTXT:
                wt.printif('... reading file ' + f)
                AllF = AllF.append(wt.get_data(f,p))
                # remove hourly text files at the end of the day
                if len(InTXT)==24 and cl.SWITCH_CLEANUP and cl.SWITCH_OUTNC:
                    os.remove(f)
        else:
            wt.printif( '... Please check SWITCH_POOL in config file!' )
        
    INPUTTIME = wt.timer(STARTTIME)

    # change scan IDs of LOS to composite VAD
    AllF = wt.change_scan_IDs(AllF)

    # export content of data frame to netcdf
    if p<>'cnr' and cl.SWITCH_OUTNC:
        wt.printif('... export nc')
        wt.export_to_netcdf(AllF,p,sDate,'')
        EXPORTTIME = wt.timer(STARTTIME)

    # plot time series (vertical line-of-sight only, 24h)
    if p<>'spectra':
        if ('LOS90' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                and cl.SWITCH_PLOT:
            if p=='wind':
                wp.plot_ts(AllF,'cnr',sDate,['dummy'])
                TSTIME = wt.timer(STARTTIME)
            wp.plot_ts(AllF,p,sDate,['dummy'])
            TSTIME = wt.timer(STARTTIME)
    
        if p<>'dbs':
            # plot low level scans (polar)
            if ('LOW' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                    and cl.SWITCH_PLOT:
                wp.plot_low_scan(AllF, p, sDate)
                if p=='wind':
                    wp.plot_low_scan(AllF,'cnr',sDate)
                LOWTIME = wt.timer(STARTTIME)

            # plot line-of-sight scans (scan duration)
            if ('LOS' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                    and cl.SWITCH_PLOT:
                wp.plot_los(AllF, p, sDate)
                if p=='wind':
                    wp.plot_los(AllF,'cnr',sDate)
                LOSTIME = wt.timer(STARTTIME)

            # calculate horizontal wind speed and direction
            if 'VAD' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE:
                if p=='wind':
                    wt.printif('... fitting radial wind')
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
            wt.printif('.. finished ' + p)
        else:
            print( '.. no output selected. Abort execution! \
                    Please adjust run options in config_lidar.py.' )

