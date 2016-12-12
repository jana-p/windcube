import datetime as dt
import matplotlib
import platform
if platform.system=='Windows':
    matplotlib.rcParams['backend'] = "Qt4Agg"
else:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.dates as dates

import numpy as np
import pandas as pd
import seaborn as sns
import pdb

# contains windcube constants
import config_lidar as cl
# contains windcube functions
import windcube_io as wio



# general plot settings
sns.set(font_scale=1.3)
sns.set_style("white")


# calculates pivot of pandas data frame, returns also axes limits 
# and color bar properties
def prepare_plotting(dfplot, sProp, pp):
    alpha = 1.0
    if pp[0]=='dummy' or pp[0]=='los':
        # for 'raw' data
        t=[x[0] for x in dfplot.index]
        if pp[0]=='dummy':
            r=[x[1] * np.sin( np.radians( dfplot.ele[0] ) ) for x in dfplot.index]
        elif pp[0]=='los':
            r=[x[1] for x in dfplot.index]
        CM='jet'
        clim1, clim2, z, CBlabel = get_lims(dfplot, sProp)
    elif pp[0]=='low_scan':
        t=np.radians(dfplot.azi)
        # distance from Mace Head (at ground)
        r=[x[1] * np.cos( np.radians( dfplot.ele[0] ) ) for x in dfplot.index]
        clim1, clim2, z, CBlabel = get_lims(dfplot, sProp)
        if sProp=='wind':
            clim1 = clim1 * 10.0
            clim2 = clim2 * 10.0
        CM='jet'
    else:
        # for wind components
        z=dfplot[pp[0]][dfplot['confidence_index']>=50]
        t=[x[0] for x in z.index]
        r=[x[1] for x in z.index]
        CBlabel=pp[1]
        wix=cl.VarDict['VAD']['cols'].index(pp[0])
        clim1=cl.VarDict['VAD']['lims'][wix][0]
        clim2=cl.VarDict['VAD']['lims'][wix][1]
        if pp[0]=='w':
            CM='coolwarm'
        elif pp[0]=='wdir':
            alpha = 0.8
            CM='hsv'
        else:
            CM='BuPu'

    # bring data from 1 dimension to a grid (2D)
    bpivot=pd.pivot(t,r,z)

    return bpivot, t, r, z, clim1, clim2, CBlabel, CM, alpha


def drop_columns(df, sProp, p):
    if p=='dummy' or p=='los' or p=='low_scan':
        # for 'raw' data
        df = df.loc[:,['ele', 'azi', cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']], 'confidence_index']]
    else:
        # for wind components
        df = df.loc[:,['ele', 'azi', p, 'confidence_index']]

    return df


def get_lims(dfplot, sProp):
    # get beta data on a logarithmic scale for plotting
    if sProp=='beta':
        z = np.log10(dfplot[cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]])
        CBlabel = cl.VarDict[sProp]['cols'][ cl.VarDict[sProp]['N'] ] \
                + ' / log10(' \
                + cl.VarDict[sProp]['units'][ cl.VarDict[sProp]['N'] ] + ')'
    else:
        z = dfplot[cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]]
        CBlabel = cl.VarDict[sProp]['cols'][ cl.VarDict[sProp]['N'] ] \
                + ' / ' \
                + cl.VarDict[sProp]['units'][ cl.VarDict[sProp]['N'] ]

    # change scale for vertical wind
    if sProp=='wind':
        clim1 = cl.VarDict[sProp]['lims'][ cl.VarDict[sProp]['N'] ][0] / 10
        clim2 = cl.VarDict[sProp]['lims'][ cl.VarDict[sProp]['N'] ][1] / 10
    else:
        clim1 = cl.VarDict[sProp]['lims'][ cl.VarDict[sProp]['N'] ][0]
        clim2 = cl.VarDict[sProp]['lims'][ cl.VarDict[sProp]['N'] ][1]

    return clim1, clim2, z, CBlabel


# plot time series
def plot_ts(AllB,sProp,sDate,plotprop):
    wio.printif('... plot ts of ' + sProp + ', ' + plotprop[0])
    # select only vertical line of sight (elevation >= 89.5)
    if plotprop[0]=='dummy':
        b1 = AllB[AllB.ele>=89.5]
        # get previously pickled data
        idlist = b1.scan_ID.unique().tolist()
        b1, dummy = wio.get_pickle(sProp, idlist, b1)
        b1 = b1[~b1.index.duplicated(keep='last')]
        b1 = b1[b1.ele>=89.5]
        # drop unnecessary columns
        b1 = drop_columns(b1, sProp, plotprop[0])
        # reduce time resolution
        b1 = b1.unstack(level='range').resample(cl.xres).mean().stack(level='range')
        name = cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']]
        title = cl.VarDict[sProp]['longs'][cl.VarDict[sProp]['N']] \
                + ' (elevation >= 89.5), on ' + sDate
    else:
        b1 = AllB
        numlim = 50.0
        # get previously pickled data
        if 'scan_ID' in b1:
            idlist = b1.scan_ID.unique().tolist()
            b1, dummy = wio.get_pickle(sProp, idlist, b1)
        elif len(plotprop)==4:
            b1, dummy = wio.get_pickle(sProp, ['VAD',plotprop[2]], b1)
        if plotprop[0]=='los':
            name = cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] \
                    + '_elev' + plotprop[2] + '_az' + plotprop[3] + '_scan'\
                    + plotprop[4]
            title = cl.VarDict[sProp]['longs'][cl.VarDict[sProp]['N']] + ' ('\
                    + plotprop[2] + ' degrees elevation), on ' + sDate
        else:
            name = plotprop[0] + '_' + plotprop[2]
            title = sProp + ' ' + ' (' + plotprop[0] + ', ' + plotprop[2] \
                    + ' degrees elevation), on ' + sDate

    # discard background (low confidence index)
    if 'confidence_index' in AllB and cl.SWITCH_REMOVE_BG:
        b = b1[b1.confidence_index>=50]
        numlim = 50.0
    else:
        b = b1
        numlim = 20.0

    # separate time and range index arrays
    b.drop_duplicates(inplace=True)
    bpivot, t, r, z, clim1, clim2, CBlabel, CM, alpha = prepare_plotting(b, sProp, plotprop)

    if cl.SWITCH_ZOOM and (sProp=='cnr'):
        limdiff = clim2 - clim1
        clim2 = clim1 + limdiff/10.0
        clim1 = clim1 - limdiff/10.0
        wio.printif('.... zooming in')
        wio.printif([clim1, clim2])
            
    # plotting
    plt.figure(figsize=(10, 5))
    plotarr = np.ma.masked_invalid(bpivot.values.T)
    cp =  plt.pcolormesh(bpivot.index.values, bpivot.columns.values, 
            plotarr, cmap=CM, edgecolors='none',
            vmin=clim1, vmax=clim2, alpha=alpha
            )
    cp.cmap.set_bad('white')

    if cl.SWITCH_ZOOM and (sProp=='cnr'):
        cp.cmap.set_over('grey')

    if plotprop[0]=='wspeed':
        cp.cmap.set_over('indigo')
    elif plotprop[0]=='w':
        cp.cmap.set_under('navy')
        cp.cmap.set_over('brown')
    cp.set_clim(clim1, clim2)
    cb = plt.colorbar(cp, extend='both')
    cb.set_label(CBlabel)
    # set axes limits and format
    # times
    h1 = int( cl.xlim[0][0:2] )
    m1 = int( cl.xlim[0][2:4] )
    s1 = int( cl.xlim[0][4:6] )
    h2 = int( cl.xlim[1][0:2] )
    m2 = int( cl.xlim[1][2:4] )
    s2 = int( cl.xlim[1][4:6] )
    StartTime = bpivot.index[0].replace(hour=h1).replace(minute=m1).replace(second=s1)
    EndTime = bpivot.index[0].replace(hour=h2).replace(minute=m2).replace(second=s2)
    plt.xlim( [StartTime, EndTime] )
    ax=plt.gca()
    if m2>0:
        hend = h2+1
    else:
        hend = h2
    if hend-h1>=6:
        ax.xaxis.set_major_locator(dates.HourLocator(byhour=range(h1,hend,(hend-h1)/6)))
    else:
        ax.xaxis.set_major_locator(dates.AutoDateLocator())
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    plt.xlabel('time / UTC')
    # range
    if plotprop[0]=='dummy':
        plt.ylim( cl.TSylim )
    elif plotprop[0]=='los' and plotprop[2]<>'90':
        if int(plotprop[4]) in cl.LOSzoom:
            plt.ylim(cl.LOSzoom[ int(plotprop[4]) ])
        else:
            plt.ylim([0,bpivot.columns[-1] * np.cos( np.radians( b.ele[0] ) )])
    elif plotprop[0]=='los' and plotprop[2]=='90':
        plt.ylim( cl.TSylim )
    else:
        plt.ylim( cl.VADylim )
    plt.ylabel('altitude agl / m')
    plt.title(title)
    plt.tight_layout()
    plt.grid(b=False)
    pdb.set_trace()
    # save and close plot
    plt.savefig(cl.figOUT + name + '_latest.png', dpi=150)
    plt.close()


# plots low level scan on polar grid
def plot_low_scan(toplot, sProp, sDate):
    if len(toplot.scan_ID)>0:
        # time difference of 59 seconds (in ns)
        newscan=np.where(np.diff(toplot.index.get_level_values('time')).astype(float) > 59000000000)
        n=0
        for s in newscan[0]:
            fig=plt.figure(figsize=(6, 5))
            # plot horizontal scan from n to s
            thisscan=toplot[n:s]
            # discard background (low confidence index)
            if 'confidence_index' in thisscan and cl.SWITCH_REMOVE_BG:
                thisscan = thisscan[thisscan.confidence_index>=50]#.status==1]#
            sTitle = 'scan on ' + thisscan.index[0][0].strftime('%Y/%m/%d')\
                    + ' from ' + thisscan.index[0][0].strftime('%H:%M:%S')\
                    + ' to ' + thisscan.index[-1][0].strftime('%H:%M:%S')
            bpivot, a, r, z, clim1, clim2, CBlabel, CM, alpha = prepare_plotting(thisscan, sProp, ['low_scan'])

            bpivotsmooth=np.zeros(np.shape(bpivot))
            window_size=5
            window = np.ones(int(window_size))/float(window_size)

            # plotting
            ax = plt.subplot(111, polar=True)
            plt.title( sTitle )
            plotarr = np.ma.masked_invalid(bpivot.values.T)
            cp =  plt.pcolormesh(bpivot.index.values, bpivot.columns.values,
                    plotarr, cmap=CM, edgecolors='none',
                    vmin=clim1, vmax=clim2, alpha=alpha
                    )
            cb = plt.colorbar(cp)
            cb.set_label(CBlabel)
            ax.set_ylim(cl.LOWylim)
            ax.set_theta_zero_location('N')
            ax.set_theta_direction(-1)
            plt.tight_layout()
            plt.grid(b=True, which='both')
            # save plot
            plt.savefig(cl.figOUT + sDate + '_' \
                    + cl.VarDict[sProp]['cols'][cl.VarDict[sProp]['N']] \
                    + '_elev_' + str( int( round(toplot.ele[0]) ) ) + '_' \
                    + str(n) + '_low_scan.png', dpi=150)
            plt.close()
            n=s


# plots LOS as time series
def plot_los(toplot, sProp, sDate):
    if len(toplot.scan_ID)>0:
        elestr = str( int( round( toplot['ele'][0] ) ) )
        azstr = str( int( round( toplot['azi'][0] ) ) )
        scanstr = str( int( round( toplot['scan_ID'][0] ) ) )
        plot_ts(toplot,sProp,sDate,['los', '', elestr, azstr, scanstr])


# plot correlations
def plot_correlation(df, p, sDate, xName, yName, sTitle, titleadd, dims):
    # set nan all outliers (out of plotting domain given by "dims")
    df[xName][ (df[xName] < dims[0]) | (df[xName] > dims[1]) ] = np.nan
    df[yName][ (df[yName] < dims[0]) | (df[yName] > dims[1]) ] = np.nan
    # plotting
    plt.figure(figsize=(7, 6))
    ax = plt.gca()
    df.plot( x=xName, y=yName, kind='hexbin', gridsize=30, ax=ax, alpha=0.8, \
            title=sTitle + titleadd)
    plt.xlabel( xName )
    plt.ylabel( yName )
    plt.xlim( dims )
    plt.ylim( dims )
    ax.plot(ax.get_xlim(), ax.get_ylim(), alpha=0.9)  # 1:1 line
    plt.tight_layout()
    plt.savefig(cl.figOUT + sDate + '_' + xName + '_' + yName + '_' \
            + titleadd + '_scatter.png', dpi=150)
    plt.close()

