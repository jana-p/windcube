import datetime as dt
import os
import pandas as pd
from glob import glob
import pdb

# contains all windcube constants:
import config_lidar as cl

# contain windcube functions:
import windcube_io as wio
import windcube_tools as wt

if cl.SWITCH_TIMER:
    import time

#import warnings
#warnings.simplefilter("ignore", FutureWarning)
#warnings.simplefilter("error")



def main(sDate,p):
    STARTTIME = time.time()

    # find files to read
    fend = cl.VarDict[p]['cols'][cl.VarDict[p]['N']]
#   InNC = sorted(glob(cl.OutPath + sDate + '_' + fend + '*.nc'))
    InTXT = sorted(glob(cl.txtInput + cl.VarDict[p]['fend'] + '.' + cl.ending))
    InPKL = sorted(glob(cl.ncInput + fend + '*.pkl'))
    for pkl in InPKL:
        os.remove(pkl)

    # read input depending on settings in config file
#   if InNC and cl.SWITCH_INPUT=='netcdf':
#       # 1) read netcdf files in list
#       wio.printif('... netcdf file found')
#       for ncf in InNC:
#           wio.printif('... reading file ' + ncf)
#           df = wio.open_existing_nc(ncf)
#           # rename variables for consistent input
#           if 'dv' in df:
#               df.rename(columns={'dv':'radial_wind_speed'}, inplace=True)
#           df = df.unstack().stack().sort_index().drop_duplicates()
#           df = df[~df.index.duplicated(keep='last')]
#           df = wt.change_scan_IDs(df)
    if cl.SWITCH_INPUT=='text':
        # 2) read all available txt files
        if InTXT:
            wio.printif('... loop over ascii files')
            # run loop if number of pools is 0 (not parallel)
            dflist = []
            for f in InTXT:
                wio.printif('... reading file ' + f)
                df = wio.get_data(f,p)
                dflist.append(df)
                if cl.SWITCH_CLEANUP:
                    wio.printif("... CLEANUP: deleting today's text files")
                    os.remove(f)
            df = pd.concat(dflist)
            df = wt.change_scan_IDs(df)
        else:
            wio.printif( '... no new text file' )
        
    INPUTTIME = wt.timer(STARTTIME)

    # export content of data frame to netcdf and pickle
    if cl.SWITCH_OUTNC and cl.SWITCH_INPUT<>'netcdf':
        wio.printif('... export nc')
        wio.export_to_netcdf(df,p,sDate,'')
        EXPORTTIME = wt.timer(STARTTIME)

    # remove yesterday's pickle and ascii files
    if cl.SWITCH_CLEANUP and cl.SWITCH_OUTNC:
        wio.printif("... CLEANUP: deleting yesterday's pickle and text files")
        wio.clean_up(sDate, fend, p)
    elif cl.SWITCH_CLEANUP and not cl.SWITCH_OUTNC:
        wio.printif('.. No cleaning up, because netcdf output is deactivated!')


if __name__=="__main__":
    # creating output path if not existing
    if os.path.exists(cl.OutPath):
        if not os.path.exists(cl.ncOUT):
            os.makedirs(cl.ncOUT)
    else:
        print("Warning: Output path doesn't exist!")

    # running "main" function, starting processing
    p = 'spectra'
    if cl.SWITCH_PLOT or cl.SWITCH_OUTNC:
        main(cl.sDate,p)
        wio.printif('.. finished ' + p)
    else:
        print( '.. no output selected. Abort execution! \
                Please adjust run options in config_lidar.py.' )

