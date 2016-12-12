import datetime as dt
import time
import os

import numpy as np
import pandas as pd
import pdb
from glob import glob
import xarray

# contains windcube constants and custom settings
import config_lidar as cl
import windcube_plotting as wp
import windcube_tools as wt



# read text files from MySQL data base, returns pandas data frame
def get_data(file_path, sProp):
    # dsb text file have a different time stamp
    if sProp=='dbs':
        dparse = lambda d: dt.datetime.strptime(d, '%Y-%m-%d %H:%M:%S')
    else:
        dparse = lambda d: dt.datetime.strptime(d, '%Y-%m-%d %H:%M:%S.%f')
    
    outdf = pd.read_csv(file_path, sep=cl.sep,          # read file from csv
            header=None, 
            skiprows=cl.skip,
            engine='python',
            names=cl.VarDict[sProp]['cols'],
            parse_dates=[0],                            # feed the first column to the parser
            date_parser=dparse, 
            squeeze=True                                # convert to `Series` object because we only have one column
            )

    outdf.drop_duplicates( inplace=True )
    # convert azimuth range from -180 - 180 degrees to 0 - 360 degrees
    outdf.loc[outdf.azi < 0, 'azi'] = outdf.loc[outdf.azi < 0, 'azi'] + 360
    # convert radial wind 'positive towards lidar' to radial wind 'positive away from lidar
    if 'radial_wind_speed' in outdf:
        outdf.radial_wind_speed = outdf.radial_wind_speed * (-1.0)

    outdf['range'] = outdf['range'].astype( float )     # change date type of range from integer to float
    outdf = outdf.set_index(['time', 'range'])          # index are time and range
    
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


# opens existing netcdf using netCDF4 and returns pandas data frame
def open_with_netcdf4(InFile):
    from netCDF4 import Dataset
    # read netcdf files
    invars = Dataset(InFile, 'r')
    # put array data to dict
    indict = {}
    namelist = []
    rangeshape = np.shape(invars.variables['range'][:])
    for var in invars.variables:
        inarr = invars.variables[var][:]
        if len(np.shape(inarr))==2:
            inarr = inarr.reshape( np.shape(inarr)[0]*np.shape(inarr)[1], )
        elif len(np.shape(inarr))==1:
            inarr = np.tile( inarr, rangeshape[1] )
        indict[var] = inarr 
        namelist.append(var)
    # create pandas data frame
    df = pd.DataFrame(columns=namelist)
    df.from_dict( indict )

    return df


# find and read appropriate pickle file
def get_pickle(p, slist, df):
    pkllist = []
    if p=='cnr':
        p = 'wind'
    if slist[0]=='VAD':
        pkln = cl.ncInput + 'VAD'
        slist = slist[1:]
    else:
        pkln = cl.ncInput + cl.VarDict[p]['cols'][cl.VarDict[p]['N']]
    for s in slist:
        pklname = pkln + '_' + str(int(s)) + '.pkl'
        InPKL = sorted(glob(pklname))
        if InPKL:
            dfp = pd.read_pickle(InPKL[0]).reset_index()
            dfp.time = pd.to_datetime(dfp.time, unit='s')
            dfp.set_index(['time','range'], inplace=True)
            if 'scan_ID' in df:
                df = pd.concat([dfp,df.loc[df.scan_ID==s]])
            else:
                df = pd.concat([dfp,df])
            df = df.drop_duplicates()
            pkllist.append(InPKL[0])

    return df, pkllist


# write all scans to individual pickle files
def write_pickle(df, pkllist, slist):
    df.drop_duplicates(inplace=True)
    for i in range(0,len(pkllist)):
        if 'scan_ID' in df and len(pkllist) == len(slist):
            df.loc[df.scan_ID==slist[i]].to_pickle(pkllist[i])
        else:
            df.to_pickle(pkllist[i])


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


# prepares pandas data frame for export to netcdf file
def export_to_netcdf(df,sProp,sDate,nameadd):
    if nameadd == '':
        sID = df.scan_ID.unique()
        for s in sID:
            dfsID = df[df['scan_ID']==s].copy()
            dfsID, pkllist = get_pickle(sProp, [s], dfsID)
            pklname = cl.ncInput + cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] + '_' + str(s) + '.pkl'
            if not pkllist:
                dfsID.to_pickle(pklname)
            else:
                dfsID.to_pickle(pkllist[0])

            # put spectra data in 3 dimensions (time, range, frequency)
            if sProp=='spectra':
                specS = dfsID.spectra.str.split(',', expand=True).astype(float).stack()
                specdf = specS.to_frame()
                # add frequency index
                specdf.reset_index(inplace=True)
                specdf.rename( columns={ 'level_2' : 'frequency_bins' }, inplace=True )
                specdf['frequency_bins'] = specdf.frequency_bins.astype(float)
                specdf.set_index(['time','range','frequency_bins'], inplace = True)
                specdf.columns = ['spectra']
            else:
                specdf = 'dummy'

            create_xarray_dataset(dfsID, nameadd, s, sProp, sDate, specdf)
#           export_with_netcdf4(dfsID, nameadd, s, sProp, sDate, specdf)
    elif 'VAD' in nameadd:
        specdf = 'dummy'
#       export_with_netcdf4(df, nameadd, 'VAD', sProp, sDate, specdf)
        create_xarray_dataset(df, nameadd, 'VAD', sProp, sDate, specdf)


# change time index to seconds since 1970 for storing in netcdf
def change_time_index(df,inds):
    df.reset_index(inplace=True)
    t = df.time.astype(np.int64) / 10**9
    df['time'] = t
    df.set_index(inds, inplace = True)

    return df


# exports xarray data set to netcdf file, including global attributes, long names and units
def create_xarray_dataset(df, nameadd, s, sProp, sDate, df3D):
    printif('.... convert from df to xarray ds')
    # change time index to seconds since 1970 for storing in netcdf
    df = change_time_index(df, ['time','range'])
    xData = xarray.Dataset.from_dataframe(df)
    xData1Ddict = {}
    # add variable attributes
    if 'VAD' in nameadd:
        pname='VAD'
    else:
        pname=sProp
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
        xData1D = xarray.Dataset(xData1Ddict).drop('range')
        xOut = xData.merge(xData1D)
        if sProp=='spectra':
            df3D = change_time_index(df3D, ['time','range','frequency_bins'])
            fbins = df3D.index.get_level_values(level='frequency_bins')
            xData3D = xarray.Dataset.from_dataframe(df3D)
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
    xOut.to_netcdf(path=cl.OutPath + sDate + '_' + nameadd + '.nc',
            mode='w', engine='netcdf4')#, format='NETCDF4')
    xOut.close()


# clean up old pickle and ascii files
def clean_up(sDate, fend, p):
    dtDate = dt.datetime.strptime( sDate, '%Y%m%d' ) - dt.timedelta(days=1)
    yesterday = dt.datetime.strftime( dtDate, '%Y%m%d')
    oldfilename = os.path.split(cl.ncInput)[0] + os.sep + yesterday
    oldpkl = sorted(glob( oldfilename + '_' + fend + '*.pkl' ))
    if oldpkl and cl.SWITCH_CLEANUP:
        wio.printif('... remove old pickle')
        [os.remove(f) for f in oldpkl]
    oldtxt = sorted(glob( oldfilename + '*' + cl.VarDict[p]['fend'] + '.txt' ))
    if oldtxt and cl.SWITCH_CLEANUP:
        wio.printif('... remove old ascii')
        [os.remove(f) for f in oldtxt]


# prints message if output option is set in config file
def printif(message):
    if cl.SWITCH_OUTPUT:
        print(message)

