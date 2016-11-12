import datetime as dt
import time
import os

import numpy as np
import pandas as pd
import pdb
import xray

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
    printif('.... convert from df to xray ds')
    if nameadd == '':
        sID = df.scan_ID.unique()
        for s in sID:
            dfsID = df[df['scan_ID']==s]
            # put spectra data in 3 dimensions (time, range, frequency)
            if sProp=='spectra':
                specS = dfsID.spectra.str.split(',', expand=True).astype(float).stack()
                specdf = specS.to_frame()
                # add frequency index
                specdf.reset_index(inplace=True)
                specdf.rename( columns={ 'level_2' : 'frequency_bins' }, inplace=True )
                specdf.set_index(['time','range','frequency_bins'], inplace = True)
                fbix = specdf.unstack().columns.labels[1]
                vals = specdf.unstack().unstack().values
                vals = vals.reshape( np.shape(vals)[0], np.shape(vals)[1]/len(fbix), len(fbix) )
            else:
                fbix = 'dummy'
                vals = 'dummy'
            create_xray_dataset(dfsID, nameadd, s, sProp, sDate, fbix, vals)
    elif 'VAD' in nameadd:
        fbix = 'dummy'
        vals = 'dummy'
        create_xray_dataset(df, nameadd, 'VAD', sProp, sDate, fbix, vals)


# exports xray data set to netcdf file, including global attributes, long names and units
def create_xray_dataset(df, nameadd, s, sProp, sDate, fbix, vals):
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


# prints message if output option is set in config file
def printif(message):
    if cl.SWITCH_OUTPUT:
        print(message)

