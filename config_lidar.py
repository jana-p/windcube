# all constants necessary to read, convert and plot the MySQL output txt files
# ScanID list, netcdf variable names, units and long names, netcdf attributes
#
##### DO NOT CHANGE #####

import os
import numpy as np

#########################

# software version
version = 1.3
vstr    = '(post STSM)'


# =============================================================================
# DATE
# un-comment the following three lines to use for near real time operation
#import datetime as dt
#td=dt.datetime.utcnow()-dt.timedelta(hours=1)
#sDate=td.strftime("%Y%m%d")
# remove following line to use for near real time operation
sDate='20160408'


# DATA PATH of input files and output netcdf files
#DataPath="/home/lidar/DATA/WindCube/"
#DataPath="//10.5.4.177/mh/WindCube/PROC/2015/"
DataPath = "C:\\Users\\JANA\\Documents\\NUIG-work\\DATA\\NUIGdata\\WindCube\\20160408\\"
#DataPath="C:\\Users\\JANA\\Documents\\NUIG-work\\MaceHead\\Instruments\\WindLidar\\data_examples\\problem\\20150623\\raw\\"
# RELATIVE DATA INPUT PATH and names using sDate
import os
ncInput = DataPath + sDate[0:4] + os.sep + sDate + '_'
txtInput = DataPath + sDate[0:4] + os.sep + sDate + '-*'


# OUTPUT PATH for figures and files
#OutPath="/home/lidar/DATA/WindCube/"
OutPath = "C:\\Users\\JANA\\Documents\\NUIG-work\\DATA\\NUIGdata\\WindCube\\20160408\\"
# RELATIVE OUTPUT PATH
ncOUT = OutPath + sDate[0:4] + os.sep
figOUT = ncOUT


# =============================================================================
### RUN OPTIONS ###

# include in list, which input text files to use ('wind' includes wind and CNR data, 'beta' includes relative backscatter, 'dbs' is for testing only)
proplist = ['wind']#,'beta','wind','dbs','spectra']

# switches
SWITCH_REMOVE_BG = True     # remove background from plot (True), or plot background (False)
SWITCH_ZOOM      = False    # zoom in to background noise (change color bar limits, only for CNR) (True), or uses limits given in VarDict (False)
SWITCH_PLOT      = True     # plot results (True)
SWITCH_OUTNC     = True     # plot results (True)
SWITCH_INNC      = False    # uses existing netcdf files if in data path (True, faulty!), or uses all text files in data path as input (False), or appends latest text file in data path to existing netcdf file in data path and removes this text file ('append', also faulty!)
SWITCH_OUTPUT    = True     # prints status messages on screen if run from command line (True)
SWITCH_TIMER     = True     # times the main processes while running the script, prints time elapsed since start of script if output is activated (True)
SWITCH_HDCP2     = False    # prepares two output files in HDCP2 format (level 1: radial wind and beta, level 2: wind components from VAD scans) (True)
SWITCH_MODE      = ['VAD','LOS90']    # calculates/plots only certain scan types ('VAD', 'LOW', 'LOS', 'LOS90'), or all scan types ('all')


# =============================================================================
### INSTRUMENT SPECIFICATIONS ###

# Scan ID definition
# add lists of all IDs that apply
ScanID={}
# any PPI scan:
ScanID['PPI']   = [20, 24, 34, 37, 38, 48, 49, 77]
# any RHI scan:
ScanID['RHI']   = [26, 32]
# any line-of-sight measurements:
ScanID['LOS']   = [20, 30, 41, 42, 44, 45, 50]
# 89.99 degrees elevation, 0 degrees azimuth:
ScanID['LOS90'] = [30, 50]
# vertical staring, for PBL detection:
ScanID['PBL']   = [20]
# beta calibration scan (low elevation, homogeneous boundary layer):
ScanID['CAL']   = [34, 38]
# VAD scan (full PPI, 0 to 360 degrees azimuth, for wind fit):
ScanID['VAD']   = [37, 48, 49, 77, 110, 111]
# any low level scans (low elevation):
ScanID['LOW']   = [20, 34, 38, 49]
# any composite of line-of-site measurements for VAD:
ScanID['COM']   = [110, 111] # 94, 

# dictionary of VAD composites
CompDict = {#94:[81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92], 
        111:[81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92], 
        110:[97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 108, 109]
        }


# =============================================================================
### OUTPUT SPECIFICATIONS ###

# LOS zoom (range axis)
# plots the specified range only, if given; otherwise plots the whole profile
# dictionary keys are scan IDs, values are a list of plot limits in meter [min_range, max_range]
LOSzoom={}
LOSzoom[41] = [5500, 6500]
LOSzoom[42] = [2000, 3000]
LOSzoom[44] = [5000, 5500]
LOSzoom[45] = [2000, 6500]


# attributes for netcdf creation
GloDict={"Location"       : 'Mace Head Atmospheric Research Station',  # optional
        "Institution"     : 'Centre for Climate and Air Pollution Studies (C-CAPS), National University of Ireland Galway, Ireland',  # Owner Institution or distributor of data set
        "Owner"           : 'Irish Aviation Authority (IAA)',  # optional
        "Title"           : 'Scanning Doppler lidar data',  # Short Title including Instrument and content of data set
        "Contact_person"  : 'Jana Preissler (jana.preissler@nuigalway.ie)',  # <Name>, <email>
        "Source"          : 'WindCube 200S, Leosphere',  # Instrument(s)
        "History"         : 'Data processed by windcube software package, version ' + str(version) + ' ' + vstr,  # How is the data set processed?
#       "Dependencies"    : 'external',  # just in case of higher level products: <file name> (without date) of the depending data set or "external" (for all data sets not archived in the data base)
        "Conventions"     : 'CF-1.6',
        "Author"          : 'Jana Preissler (jana.preissler@nuigalway.ie)',
        "Comments"        : '',  # Miscellaneous Information about your dataset
        "Licence"         : 'For non-commercial use only.',  #  This data is subject to the HD(CP)2 data policy to be found at www.hdcp2.eu and in the HD(CP)2 Observation Data Product standard.
        "year"            : np.int16( sDate[0:4] ),  # year of dataset
        "month"           : np.int16( sDate[4:6] ),  # month of dataset
        "day"             : np.int16( sDate[6:] )  #  dataset
        }

# general variables additional to those given in VarDict (added to all netCDF files)
GenDict={"cols"  : ("lat", "lon", "zsl", "wl"),
        "names"  : ("latitude", "longitude", "altitude", "radiation_wavelength"),
        "units"  : ("degree_north", "degree_east", "m", "m"),
        "longs"  : ("latitude of sensor", "longitude of sensor", "altitude of sensor above mean sea level", "wavelength"),
        "val"    : (53.33, -9.9, 21.0, 0.000001543),
        "ty"     : ('d', 'd', 'd', 'd')
        }



# variables for plotting and netcdf creation:
# there is one entry in the dictionary for each output parameter (level 0: spectra; level 1: wind, cnr, beta; level 2: VAD)
# fend  ... file name ending
# N     ... position of variable to plot in text file
# cols  ... column names in text file
# names ... variable names for output files (netcdf, HD(CP)2 convention, version 1.5)
# units ... units for output files (netcdf, HD(CP)2 convention, version 1.5)
# longs ... long names for output files (netcdf, HD(CP)2 convention, version 1.5)
# lims  ... axis limits for plotting
# ty    ... data type (for netcdf creation?)
VarDict={
        "wind" : {"fend"  : "_whole_radial",
                  "N"     : 7,
                  "cols"  : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "names" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "m s-1", "m s-1", "dB", "no unit", "percent", "no unit"),
                  "longs" : ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "radial wind speed", "standard deviation of the radial wind speed", "carrier-to-noise ratio", "confidence index", "mean error", "status"),
                  "lims"  : ([], [], [], [], [-180,180], [0,90], [0,8000], [-20,20], [-1,1], [-30,10], [0,100], [], [0,1]),
                  "ty"    : ('d', 'i', 'i', 'i', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'b')
            },
        "VAD"  : {"fend"  : "_whole_radial",
                  "N"     : 2,
                  "cols"  : ("time", "range", 'wspeed', 'w', 'wdir', 'confidence_index', 'rsquared', 'number_of_function_calls'),
                  "names" : ("time", "range", 'wspeed', 'w', 'wdir', 'confidence_index', 'rsquared', 'number_of_function_calls'),
                  "units" : ("seconds since 1970-01-01 00:00:00", "meter above ground level", "m s-1", "m s-1", "degree", "no unit", "no unit", "no unit"),
                  "longs" : ("time", "distance from sensor to center of each range gates along the line of sight", "wind_speed", "upward_air_velocity", "wind_from_direction", "confidence index", "R^2 of fit", "number of function calls for least spuare fit"),
                  "lims"  : ([], [0,8000], [0,25], [-2,2], [0,360], [0,100], [0,1], []),
                  "ty"    : ('d', 'd', 'd', 'd', 'd', 'd', 'd', 'b')
            },
        "cnr"  : {"fend"  : "_whole_radial",
                  "N"     : 9,
                  "cols"  : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "names" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "m s-1", "m s-1", "dB", "no unit", "percent", "no unit"),
                  "longs" : ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "radial wind speed", "standard deviation of the radial wind speed", "carrier-to-noise ratio", "confidence index", "mean error", "status"),
                  "lims"  : ([], [], [], [], [-180,180], [0,90], [0,8000], [-20,20], [-1,1], [-30,10], [0,100], [], [0,1]),
                  "ty"    : ('d', 'i', 'i', 'i', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'b')
            },
        "beta" : {"fend"  : "_radial_beta",
                  "N"     : 7,
                  "cols"  : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "beta"),
                  "names" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "beta"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "1/(m*sr)"),
                  "longs" : ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "attenuated relative backscatter"),
                  "lims"  : ([], [], [], [], [-180,180], [0,90], [0,8000], [-8,-3]),
                  "ty"    : ('d', 'i', 'i', 'i', 'd', 'd', 'd', 'd')
            },
        "dbs"  : {"fend"  : "_dbs_wind",
                  "N"     : 6,
                  "cols"  : ("time", "azi", "ele", "range", "xwind", "ywind", "zwind", "CNR", "confidence_index"),
                  "names" : ("time", "azi", "ele", "range", "u", "v", "upward_air_velocity", "CNR", "confidence_index"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "degree", "degree", "m", "m s-1", "m s-1", "m s-1", "dB", "1"),
                  "longs" : ("time", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "horizontal wind (east-west)", "horizontal wind (north-south)", "Vertical velocity", "carrier-to-noise ratio", "confidence index"),
                  "lims"  : ([], [-180,180], [0,90], [0,8000], [-20,20], [-20,20], [-2,2], [-30,10], [0,100]),
                  "ty"    : ('d', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd')
            },
        "spectra":{"fend" : "_spectra",
                   "N"    : 7,
                   "cols" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "spectra"),
                   "names": ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "spectra"),
                   "units": ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "a.u."),
                   "longs": ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "spectral density"),
                   "lims" : ([], [], [], [], [-180,180], [0,90], [0,8000], [0,5000]),
                   "ty"   : ('d', 'i', 'i', 'i', 'd', 'd', 'd', 'd')
            }
        }


# variable attributes: [name (in netcdf file), standard name, long name, unit, type, comments, dimensions, hdcp2]
AttDict={
        "time"                  : ["time", "time", "time", "seconds since 1970-01-01 00:00:00", "d", "", [], True],
        "range"                 : ["range", "", "distance from sensor to center of range gate along the line of sight", "m", "d", "", [], True],
        "azi"                   : ["azi", "sensor_azimuth_angle", "sensor azimuth due north", "degree", "d", "The reference direction is due north. The angle is measured clockwise positive, starting from north.", 1, True],
        "ele"                   : ["ele", "sensor_elevation_angle", "sensor elevation angle", "degree", "d", "Elevation angle from the horizon increasing towards zenith.", 1, True],
        "radial_wind_speed"     : ["dv", "radial_velocity_of_scatterers_away_from_instrument", "doppler velocity", "m s-1", "d", "A velocity is a vector quantity. 'Radial velocity away from instrument' means the component of the velocity of the scatterers along the line of sight of the instrument, where positive implies movement away from the instrument", 2, True],
        "dv"                    : ["dv", "radial_velocity_of_scatterers_away_from_instrument", "doppler velocity", "m s-1", "d", "A velocity is a vector quantity. 'Radial velocity away from instrument' means the component of the velocity of the scatterers along the line of sight of the instrument, where positive implies movement away from the instrument", 2, True],
        "beta"                  : ["beta", "volume_attenuated_backwards_scattering_function_in_air", "attenuated backscatter coefficient", "m-1 sr-1", "d", "determined from SNR; uncalibrated and uncorrected", 2, True],
        "att_rel_beta"          : ["beta", "volume_attenuated_backwards_scattering_function_in_air", "attenuated backscatter coefficient", "m-1 sr-1", "d", "determined from SNR; uncalibrated and uncorrected", 2, True],
        "wdir"                  : ["wdir", "wind_from_direction", "", "degree", "d", "", 2, True],
        "wspeed"                : ["wspeed", "wind_speed", "", "m s-1", "d", "", 2, True],
        "w"                     : ["w", "upward_air_velocity", "", "m s-1", "d", "A velocity is a vector quantity. 'Upward' indicates a vector component which is positive when directed upward (negative downward)", 2, True],
        "CNR"                   : ["CNR", "", "carrier-to-noise ratio", "dB", "d", "", 2, False],
        "spectra"               : ["spectra", "", "", "1", "d", "", 3, False],
        "frequency_bins"        : ["frequency_bins", "", "", "1", "d", "", [], False],
        "dev_radial_wind_speed" : ["dev_radial_wind_speed", "", "standard deviation of the radial wind speed", "m s-1", "d", "", 2, False],
        "mean_error"            : ["mean_error", "", "mean error", "percent", "d", "", 2, False],
        "confidence_index"      : ["confidence_index", "", "confidence index", "1", "d", "", 2, False],
        "config_ID"             : ["config_ID", "", "", "1", "i", "", 1, False],
        "scan_ID"               : ["scan_ID", "", "", "1", "i", "", 1, False],
        "LOS_ID"                : ["LOS_ID", "", "", "1", "i", "", 1, False],
        }
