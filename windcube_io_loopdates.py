import datetime as dt
import time
import os

import numpy as np
import pandas as pd
import pdb
from glob import glob
import xarray
from netCDF4 import Dataset

# contains windcube constants and custom settings
import config_lidar_loop as cl
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
    xds = xarray.open_dataset(InFile, decode_times=True ).transpose()
    # convert to pandas data frame
    dfnc = xds.to_dataframe()
    # swap order of indices (first time, then range)
    dfswap = dfnc.swaplevel('time','range')
    xds.close()

    return dfswap


# opens existing netcdf using netCDF4 and returns pandas data frame
def open_with_netcdf4(InFile, sDate):
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
    df = df.from_dict( indict )
    hour = df['time'].astype(int)
    minute = ((df['time'] - hour)*60).astype(int)
    second = ((df['time'] - hour - minute/60)*60*60).astype(int)
    df['time'] = [pd.Timestamp(dt.datetime(int(sDate[0:4]), int(sDate[4:6]), int(sDate[6:8]), hour.values[t], minute.values[t], second.values[t])) for t in range(0,len(hour.values))]
    # change date type of range from integer to float
    df['range'] = df['range'].astype( float )
    df = df.set_index(['time', 'range'])
    df.drop_duplicates( inplace=True )

    # convert azimuth range from -180 - 180 degrees to 0 - 360 degrees
    df.loc[df.azimuth < 0, 'azimuth'] = df.loc[df.azimuth < 0, 'azimuth'] + 360
    # convert radial wind 'positive towards lidar' to radial wind 'positive away from lidar
    if 'radial_wind_speed' in df:
        df.radial_wind_speed = df.radial_wind_speed * (-1.0)

    return df


# find and read appropriate pickle file
def get_pickle(p, slist, df, sDate):
    pkllist = []
    if p=='cnr':
        p = 'wind'
    if slist[0]=='VAD':
        pkln = cl.ncInput + sDate + 'VAD'
        slist = slist[1:]
    else:
        pkln = cl.ncInput + sDate + cl.VarDict[p]['cols'][cl.VarDict[p]['N']]
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
        if vname<>vnamenew:
            df.rename( name_dict={ vname : vnamenew }, inplace=True )
        try:
            if cl.AttDict[vname][6]==1:
                df1D[vnamenew] = df[vnamenew].mean(axis=1,skipna=True)#df[vnamenew][:,0]
                df1D[vnamenew].attrs = varatts
                df = df.drop( [vnamenew] )
            else:
                df[vnamenew].attrs = varatts
        except KeyError:
            for cols in df.keys(): 
                if vname in cols[0:3]:
                    if cl.AttDict[vname][6]==1:
                        df1D[vnamenew] = df[cols].mean(axis=1,skipna=True)#df[cols][:,0]
                        df1D[vnamenew].attrs = varatts
                        df = df.drop( [cols] )
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
            dfsID, pkllist = get_pickle(sProp, [s], dfsID, sDate)
            pklname = cl.ncInput + sDate + cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] + '_' + str(s) + '.pkl'
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
    if sProp=='spectra':
        hour1 = df.index.get_level_values(level=0)[0].hour
        hour2 = df.index.get_level_values(level=0)[-1].hour
    df = change_time_index(df, ['time','range'])
    if sProp=='wind':
        df = df[(df.CNR!=-9999) & (df.radial_wind_speed!=9999)]
    elif sProp=='beta':
        df = df[(df.beta!=-9999)]
    df.index.drop_duplicates(keep='last')
#   df[df.index.duplicated()]
    xData = xarray.Dataset.from_dataframe(df)
    xData = xData.isel(range=np.isnan(xData.CNR).sum(axis=0)<np.shape(xData.CNR)[0]-np.shape(xData.CNR)[0]/10)
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
        try:
            xData1D = xarray.Dataset(xData1Ddict).drop('range')
        except ValueError:
            xData1D = xarray.Dataset(xData1Ddict)
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
    if 'spectra' in nameadd:
        pathname = cl.OutPath + sDate + '_' + str(hour1) + '-' + str(hour2) + '_' + nameadd + '.nc'
    else:
        pathname = cl.OutPath + sDate + '_' + nameadd + '.nc'
    xOut.to_netcdf(path=pathname,
            mode='w', engine='netcdf4')#, format='NETCDF4')
    xOut.close()


# uses netcdf4 to write data to netcdf file (instead of xarray)
def export_with_netcdf4(df, nameadd, s, sProp, sDate, df3D):
    printif('.... export with netcdf4')
    df.reset_index(inplace=True)
    if 'VAD' in nameadd:
        pname='VAD'
    else:
        pname=sProp
    # open netcdf file to write
    if sProp=='dbs':
        nameadd = 'DBS'
    elif sProp=='hdcp2':
        nameadd = 'HDCP2_' + nameadd
    else:
        if 'VAD' in nameadd:
            nameadd = nameadd[1:]
        else:
            nameadd = cl.VarDict[pname]['cols'][cl.VarDict[pname]['N']] + nameadd + '_scanID_' + str(s)
    rootgrp = Dataset(cl.OutPath + sDate + '_' + nameadd + '.nc', 'w', format="NETCDF4")
    # define dimensions
    t = df['time'].unique().astype(np.int64) / 10**9
    r = df['range'].unique()
    df['t'] = df['time']
    df['r'] = df['range']
    df.set_index(['t','r'], inplace=True)
    tdim = rootgrp.createDimension('time', None)
    rdim = rootgrp.createDimension('range', len(r))
    # add variables and attributes
    DimDict = [('dummy',), ('time',), ('time','range')]
    if sProp=='hdcp2':
        for col in df:
            ################### change this
            var = rootgrp.createVariable(col,"f8",('time',))
#           rootgrp = all_att_to_netcdf4( rootgrp, col, pname, [] )
            ###################
    else:
        for c in range( 0, len(cl.VarDict[pname]['cols']) ):
            if not cl.AttDict[cl.VarDict[pname]['cols'][c]][6]:
                if cl.VarDict[pname]['cols'][c]=='time':
                    dims = ('time',)
                    vals = t
                elif cl.VarDict[pname]['cols'][c]=='range':
                    dims = ('range',)
                    vals = r
            else:
                dims = DimDict[cl.AttDict[cl.VarDict[pname]['cols'][c]][6]]
                if cl.AttDict[cl.VarDict[pname]['cols'][c]][6]==1:
                    vals = df[cl.VarDict[pname]['cols'][c]].unstack().values[:,0]
                elif cl.AttDict[cl.VarDict[pname]['cols'][c]][6]==2:
                    vals = df[cl.VarDict[pname]['cols'][c]].unstack().values
            # define variable
            var = rootgrp.createVariable(cl.VarDict[pname]['cols'][c],cl.VarDict[pname]['ty'][c],dims)
            # add values
            var[:] = vals
            # add variable attributes
            var.units = cl.AttDict[cl.VarDict[pname]['cols'][c]][3]
            var.long_name = cl.AttDict[cl.VarDict[pname]['cols'][c]][2]
            var.standard_name = cl.AttDict[cl.VarDict[pname]['cols'][c]][1]
            var.comments = cl.AttDict[cl.VarDict[pname]['cols'][c]][5]
    # add general variables
    for gen in cl.GenDict['cols']:
        genvar = rootgrp.createVariable(gen, cl.GenDict['ty'][cl.GenDict['cols'].index(gen)])
        genvar[:] = cl.GenDict['val'][cl.GenDict['cols'].index(gen)]
        genvar.units = cl.GenDict['units'][cl.GenDict['cols'].index(gen)]
        genvar.long_name = cl.GenDict['longs'][cl.GenDict['cols'].index(gen)]
    # add global attributes
    for glo in cl.GloDict:
        rootgrp.gloatt = cl.GloDict[glo]
        rootgrp.renameAttribute('gloatt',glo)
    rootgrp.close()
    df.index.rename(['time','range'],inplace=True)
    df.drop(['time','range'], axis=1, inplace=True)


# clean up old pickle and ascii files
def clean_up(sDate, fend, p):
    # find files
    dtDate = dt.datetime.strptime( sDate, '%Y%m%d' ) - dt.timedelta(days=1)
    yesterday = dt.datetime.strftime( dtDate, '%Y%m%d')
    oldfilename = os.path.split(cl.ncInput)[0] + os.sep + yesterday
    oldpkl = sorted(glob( oldfilename + '_' + fend + '*.pkl' ))
    oldtxt = sorted(glob( oldfilename + '*' + cl.VarDict[p]['fend'] + '.txt' ))
    if 'wind' in fend:
        oldVAD = sorted(glob( oldfilename + '_VAD' + '*.pkl' ))
    else:
        oldVAD = []
    # remove pickle files
    if oldpkl and cl.SWITCH_CLEANUP:
        printif('... remove old pickle')
        [os.remove(f) for f in oldpkl]
    # remove text files
    if oldtxt and cl.SWITCH_CLEANUP:
        printif('... remove old ascii')
        [os.remove(f) for f in oldtxt]
    # remove VAD pickle files
    if oldVAD and cl.SWITCH_CLEANUP:
        printif('... remove old VAD pickle')
        [os.remove(f) for f in oldVAD]


# prints message if output option is set in config file
def printif(message):
    if cl.SWITCH_OUTPUT:
        print(message)

