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
import xray

import config_lidar as cl                               # contains all constants
import windcube_tools as wt


sns.set(font_scale=1.3)
sns.set_style("white")


# calculates pivot of pandas data frame, returns also axes limits and color bar properties
def prepare_plotting(dfplot, sProp, pp):
    if pp[0]=='dummy' or pp[0]=='los':
        # for 'raw' data
        t=[x[0] for x in dfplot.index]
        r=[x[1] for x in dfplot.index]
        CM='jet'
        clim1, clim2, z, CBlabel = get_lims(dfplot, sProp)
    elif pp[0]=='low_scan':
        t=np.radians(dfplot.azimuth)
        r=[x[1] * np.cos( np.radians( dfplot.elevation[0] ) ) for x in dfplot.index] # distance from Mace Head (at ground)
        clim1, clim2, z, CBlabel = get_lims(dfplot, sProp)
        if sProp=='wind':
            clim1 = clim1 * 10.0
            clim2 = clim2 * 10.0
        CM='rainbow'
    else:
        # for wind components
        z=dfplot[pp[0]][dfplot['confidence_index']>=50]
        t=[x[0] for x in z.index]
        r=[x[1] for x in z.index]
        CBlabel=pp[1]
        if pp[0]=='vertical':
            CM='coolwarm'
            clim1=-2
            clim2=2
        elif pp[0]=='direction':
            CM='rainbow'
            clim1=0
            clim2=360
        else:
            CM='BuPu'
            clim1=0
            clim2=25

    # bring data from 1 dimension to a grid (2D)
    bpivot=pd.pivot(t,r,z)

    return bpivot, t, r, clim1, clim2, CBlabel, CM


def get_lims(dfplot, sProp):
    # get beta data on a logarithmic scale for plotting
    if sProp=='beta':
        z=np.log10(dfplot[cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]])
        CBlabel=cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] + ' / log10(' + cl.VarDict[sProp]['units'][cl.VarDict[sProp]['N']] + ')'
    else:
        z=dfplot[cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]]
        CBlabel=cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] + ' / ' + cl.VarDict[sProp]['units'][cl.VarDict[sProp]['N']]

    # change scale for vertical wind
    if sProp=='wind':
        clim1=cl.VarDict[sProp]['lims'][cl.VarDict[sProp]['N']][0]/10
        clim2=cl.VarDict[sProp]['lims'][cl.VarDict[sProp]['N']][1]/10
    else:
        clim1=cl.VarDict[sProp]['lims'][cl.VarDict[sProp]['N']][0]
        clim2=cl.VarDict[sProp]['lims'][cl.VarDict[sProp]['N']][1]

    return clim1, clim2, z, CBlabel


# plot time series
def plot_ts(AllB,sProp,sDate,plotprop):
    wt.printif('... plot ts of ' + sProp + ', ' + plotprop[0])
    # select only vertical line of sight (elevation >= 89.5)
    if plotprop[0]=='dummy':
        b1 = AllB[AllB.elevation>=89.5]
        name = cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]
        title = cl.VarDict[sProp]['longs'][cl.VarDict[sProp]['N']] + ' (elevation >= 89.5), on ' + sDate
        # discard background (low confidence index)
        if 'confidence_index' in AllB and cl.SWITCH_REMOVE_BG:
            b = b1[b1.confidence_index>=30]
        else:
            b = b1
    else:
        b = AllB
        if plotprop[0]=='los':
            name = cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] + '_elev' + plotprop[2] + '_az' + plotprop[3] + '_scan' + plotprop[4]
            title = cl.VarDict[sProp]['longs'][cl.VarDict[sProp]['N']] + ' (' + plotprop[2] + ' degrees elevation), on ' + sDate
        else:
            name=plotprop[0] + '_' + plotprop[2]
            title=sProp + ' ' + plotprop[0] + ' (' + plotprop[2] + ' degrees elevation), on ' + sDate

    # separate time and range index arrays
    bpivot, t, r, clim1, clim2, CBlabel, CM = prepare_plotting(b, sProp, plotprop)
    if cl.SWITCH_ZOOM and (sProp=='cnr'):
        limdiff = clim2 - clim1
        clim2 = clim1 + limdiff/10.0
        clim1 = clim1 - limdiff/10.0
        wt.printif('.... zooming in')
        wt.printif(clim1, clim2)

    # plotting
    plt.figure(figsize=(10, 5))
    cp = plt.contourf(bpivot.index, bpivot.columns, bpivot.T, cmap=CM, 
            vmin=clim1, vmax=clim2,
            levels=np.arange(clim1, clim2, (clim2-clim1)/50.0)
            )
    cb = plt.colorbar(cp)
    cb.set_label(CBlabel)
    # set axes limits and format
    if plotprop[0]=='los':
        plt.xlim([bpivot.index[0], bpivot.index[-1]])
    else:
        plt.xlim([bpivot.index[0], bpivot.index[0]+dt.timedelta(hours=24)])
    ax=plt.gca()
    ax.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    plt.title(title)
    if plotprop[0]=='dummy':
        plt.ylim([0,15000])
        plt.ylabel('altitude agl / m')
    elif plotprop[0]=='los' and plotprop[2]<>'90':
        if int(plotprop[4]) in cl.LOSzoom:
            plt.ylim(cl.LOSzoom[ int(plotprop[4]) ])
        else:
            plt.ylim([0,bpivot.columns[-1] * np.cos( np.radians( b.elevation[0] ) )])
        plt.ylabel('range / m')
    elif plotprop[0]=='los' and plotprop[2]=='90':
        plt.ylim([0,10000])
        plt.ylabel('altitude agl / m')
    else:
        plt.ylim([0,5000])
        plt.ylabel('altitude agl / m')
    plt.xlabel('time / UTC')
    plt.tight_layout()
    plt.grid(b=False)
    # save plot
    plt.savefig(cl.OutPath + sDate[0:4] + os.sep + name + '_latest.png', dpi=150)


# plots low level scan on polar grid
def plot_low_scan(AllB, sProp, sDate):
    for LOWscan in cl.ScanID['LOW']:
        # select only VAD scans
        toplot = AllB[AllB.scan_ID==LOWscan]
        if len(toplot.scan_ID)>0:
            newscan=np.where(np.diff(toplot.index.get_level_values('time')).astype(float) > 59000000000) # time difference of 59 seconds (in ns)
            n=0
            for s in newscan[0]:
                fig=plt.figure(figsize=(6, 5))
                # plot horizontal scan from n to s
                thisscan=toplot[n:s]
        
                bpivot, a, r, clim1, clim2, CBlabel, CM = prepare_plotting(thisscan, sProp, ['low_scan'])

                bpivotsmooth=np.zeros(np.shape(bpivot))
                window_size=5
                window = np.ones(int(window_size))/float(window_size)

                # plotting
                ax = plt.subplot(111, polar=True)
                cp = plt.contourf(bpivot.index, bpivot.columns, bpivot.T, cmap=CM, 
                        vmin=clim1, vmax=clim2,
                        levels=np.arange(clim1, clim2, (clim2-clim1)/50.0)
                        )
                cb = plt.colorbar(cp)
                cb.set_label(CBlabel)
                ax.set_theta_zero_location('N')
                ax.set_theta_direction(-1)
                plt.tight_layout()
                plt.grid(b=True, which='both')
                # save plot
                plt.savefig(cl.OutPath + sDate[0:4] + os.sep + sDate + '_' + 
                        cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] + 
                        '_elev_' + str( int( round(toplot.elevation[0]) ) ) + '_' + str(n) + '_low_scan.png', dpi=150)
                n=s


# plots LOS as time series
def plot_los(AllX, sProp, sDate):
    for LOSscan in cl.ScanID['LOS']:
        # select only VAD scans
        toplot = AllX[AllX.scan_ID==LOSscan]
#       pdb.set_trace()
        if len(toplot.scan_ID)>0:
            newscan=np.where(np.diff(toplot.index.get_level_values('time')).astype(float) > 59000000000) # time difference of 59 seconds (in ns)
            n=0
            for s in newscan[0]:
                elestr = str( int( round( toplot['elevation'][0] ) ) )
                azstr = str( int( round( toplot['azimuth'][0] ) ) )
                scanstr = str( int( round( toplot['scan_ID'][0] ) ) )
                plot_ts(toplot,sProp,sDate,['los', '', elestr, azstr, scanstr])

