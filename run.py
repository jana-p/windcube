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
    InNC = sorted(glob(cl.OutPath + sDate + '_' + fend + '*.nc'))
    InTXT = sorted(glob(cl.txtInput + cl.VarDict[p]['fend'] + '.' + cl.ending))

    # read input depending on settings in config file
    if InNC and cl.SWITCH_INPUT=='netcdf':
        # 1) read netcdf files in list
        wio.printif('... netcdf file found')
        for ncf in InNC:
            wio.printif('... reading file ' + ncf)
            df = wio.open_existing_nc(ncf)
            # rename variables for consistent input
            if 'dv' in df:
                df.rename(columns={'dv':'radial_wind_speed'}, inplace=True)
            df = df.unstack().stack().sort_index().drop_duplicates()
            df = df[~df.index.duplicated(keep='last')]
            df15min = run_filewise(df,p,sDate,STARTTIME,fend)
    elif (cl.SWITCH_INPUT=='text' or cl.SWITCH_INPUT=='pickle'):
        # 4) read all available txt files
        if InTXT:
            wio.printif('... loop over ascii files')
            df15min = pd.DataFrame()
            # run loop if number of pools is 0 (not parallel)
            for f in InTXT:
                wio.printif('... reading file ' + f)
                df = wio.get_data(f,p)
                df = pd.concat([df15min, df])
                df15min = run_filewise(df,p,sDate,STARTTIME,fend)
                if cl.SWITCH_CLEANUP and cl.SWITCH_OUTNC:
                    os.remove(f)
        else:
            wio.printif( '... no new text file' )

    INPUTTIME = wt.timer(STARTTIME)


def run_filewise(AllF,p,sDate,STARTTIME,fend):
    # change scan IDs of LOS to composite VAD
    AllF = wt.change_scan_IDs(AllF)
    scanIDs = AllF['scan_ID'].unique()

    # export content of data frame to netcdf and pickle
    if cl.SWITCH_OUTNC and cl.SWITCH_INPUT<>'netcdf':
        wio.printif('... export nc')
        wio.export_to_netcdf(AllF,p,sDate,'')
        EXPORTTIME = wt.timer(STARTTIME)

    # plot time series (vertical line-of-sight only, 24h)
    if p<>'spectra':
        if ('LOS90' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                and cl.SWITCH_PLOT \
                and any([a in cl.ScanID['LOS90'] for a in scanIDs]):
            if p=='wind':
                wp.plot_ts(AllF,'cnr',sDate,['dummy'])
            wp.plot_ts(AllF,p,sDate,['dummy'])
            TSTIME = wt.timer(STARTTIME)

        if p<>'dbs':
            # plot low level scans (polar)
            if ('LOW' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                    and cl.SWITCH_PLOT \
                    and any([a in cl.ScanID['LOW'] for a in scanIDs]):
                try:
                    [wp.plot_low_scan( AllF[AllF.scan_ID==LOWscan], p, sDate ) \
                            for LOWscan in cl.ScanID['LOW']]
                    if p=='wind':
                        [wp.plot_low_scan( AllF[AllF.scan_ID==LOWscan], 'cnr', sDate ) \
                                for LOWscan in cl.ScanID['LOW']]
                except ValueError:
                    wio.printif('... ValueError in LOW ' + str(AllF.scan_ID.unique()))

                LOWTIME = wt.timer(STARTTIME)

            # plot line-of-sight scans (scan duration)
            if ('LOS' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE) \
                    and cl.SWITCH_PLOT \
                    and any([a in cl.ScanID['LOS'] for a in scanIDs]):
                [wp.plot_los( AllF[AllF.scan_ID==LOSscan], p, sDate ) \
                        for LOSscan in cl.ScanID['LOS']]
                if p=='wind':
                    [wp.plot_los( AllF[AllF.scan_ID==LOSscan], 'cnr', sDate ) \
                            for LOSscan in cl.ScanID['LOS']]
                LOSTIME = wt.timer(STARTTIME)

            # calculate horizontal wind speed and direction
            if 'VAD' in cl.SWITCH_MODE or 'all' in cl.SWITCH_MODE:
                if p=='wind' and len(scanIDs)>=1 \
                    and any([a in cl.ScanID['VAD'] for a in scanIDs]):
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

    # remove yesterday's pickle and ascii files
    if cl.SWITCH_CLEANUP and cl.SWITCH_OUTNC:
        wio.clean_up(sDate, fend, p)
    elif cl.SWITCH_CLEANUP and not cl.SWITCH_OUTNC:
        wio.printif('.. No cleaning up, because netcdf output is deactivated!')

    # return last 15 minutes of hour to avoid cut-off of VAD scans
    last15 = AllF.reset_index().set_index('time').last('15T')
    last15 = last15.reset_index().set_index(['time','range'])

    return last15



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

