import datetime as dt
import time
import os

import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
#matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.dates as dates

import numpy as np
from scipy import optimize
import pandas as pd
import seaborn as sns
import pdb

import config_lidar as cl                               # contains all constants
import windcube_tools as wt


sns.set(font_scale=1.3)
sns.set_style("white")



# compares DBS wind to VAD scan fit results
def compare_dbs(DBSdf, p, sDate):
    # get VAD data (requires existing VAD netcdf file)
    VADfile = cl.OutPath + sDate[0:4] + os.sep + sDate + '_' + cl.VarDict['VAD']['cols'][cl.VarDict['VAD']['N']] + '_VAD' + '_75.nc'
    VADdf = wt.open_existing_nc(VADfile)
    VADdfix = VADdf.reset_index()
    DBSdfix = DBSdf.reset_index()
    df = pd.merge(VADdfix, DBSdfix, how='outer')
    df = df.set_index( ['time', 'range'] )

    # obtain horizontal wind speed and direction from xwind and ywind
    df['hwind'] = np.sqrt( df.xwind**2 + df.ywind**2 )
    df['hdir'] = np.degrees( np.arctan2( df.xwind, df.xwind ) ) + 180

    # test different averaging
    df1min = df.unstack('range').resample('1T', fill_method='bfill').stack('range')
    df2min = df.unstack('range').resample('2T', fill_method='bfill').stack('range')
    df3min = df.unstack('range').resample('3T', fill_method='bfill').stack('range')
    df5min = df.unstack('range').resample('5T', fill_method='bfill').stack('range')

    # plot correlation of horizontal wind speed
    plot_correlation( df1min, p, sDate, 'speed', 'hwind', 'Horizontal wind speed, ', '1 minute', [0, 30] )
    plot_correlation( df2min, p, sDate, 'speed', 'hwind', 'Horizontal wind speed, ', '2 minutes', [0, 30] )
    plot_correlation( df3min, p, sDate, 'speed', 'hwind', 'Horizontal wind speed, ', '3 minutes', [0, 30] )
    plot_correlation( df5min, p, sDate, 'speed', 'hwind', 'Horizontal wind speed, ', '5 minutes', [0, 30] )

    # plot correlation of horizontal wind direction
    plot_correlation( df1min, p, sDate, 'direction', 'hdir', 'Horizontal wind direction, ', '1 minute', [0, 360] )
    plot_correlation( df3min, p, sDate, 'direction', 'hdir', 'Horizontal wind direction, ', '3 minutes', [0, 360] )

    # plot correlation of vertical velocity
    plot_correlation( df1min, p, sDate, 'vertical', 'zwind', 'Vertical velocity, ', '1 minute', [-4, 4] )
    plot_correlation( df3min, p, sDate, 'vertical', 'zwind', 'Vertical velocity, ', '3 minutes', [-4, 4] )


# plot correlations
def plot_correlation(df, p, sDate, xName, yName, sTitle, titleadd, dims):
    # set nan all outliers (out of plotting domain given by "dims")
    df[xName][ (df[xName] < dims[0]) | (df[xName] > dims[1]) ] = np.nan
    df[yName][ (df[yName] < dims[0]) | (df[yName] > dims[1]) ] = np.nan
    # plotting
    plt.figure(figsize=(7, 6))
    ax = plt.gca()
    df.plot( x=xName, y=yName, kind='hexbin', gridsize=30, ax=ax, alpha=0.8, title=sTitle + titleadd)
    plt.xlabel( xName )
    plt.ylabel( yName )
    plt.xlim( dims )
    plt.ylim( dims )
    ax.plot(ax.get_xlim(), ax.get_ylim(), alpha=0.9)  # 1:1 line
    plt.tight_layout()
    plt.savefig(cl.OutPath + sDate[0:4] + os.sep + sDate + '_' + xName + '_' + yName + '_' + titleadd + '_scatter.png', dpi=150)

