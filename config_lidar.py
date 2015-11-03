# all constants necessary to read, convert and plot the MySQL output txt files
# ScanID list, netcdf variable names, units and long names, netcdf attributes


# data path of input files
#DataPath="/home/lidar/DATA/WindCube/"
#DataPath="//10.5.4.177/mh/WindCube/PROC/2015/"
DataPath="C:\\Users\\JANA\\Documents\\NUIG-work\\DATA\\NUIGdata\\WindCube\\20151027\\"
#DataPath="C:\\Users\\JANA\\Documents\\NUIG-work\\MaceHead\\Instruments\\WindLidar\\data_examples\\problem\\20150623\\raw\\"

# output path
#OutPath="/home/lidar/DATA/WindCube/"
OutPath="C:\\Users\\JANA\\Documents\\NUIG-work\\DATA\\NUIGdata\\WindCube\\20151027\\out_git\\"

# date
# un-comment the following two lines to use for near real time operation
#td=datetime.utcnow()-timedelta(hours=1)
#sDate=td.strftime("%Y%m%d")
# remove following line to use for near real time operation
sDate='20151027'

# include in list, which input text files to use ('wind' includes wind and CNR data, 'beta' includes relative backscatter)
proplist=['wind']#,'beta',''wind']

#### shell scripts
ShellScriptFile="/home/lidar/SCRIPTS/SH/get_WindCube_data_mysql_hourly.SH"
ShellScriptLog="/home/lidar/SCRIPTS/SH/get_WindCube_data_mysql_hourly.LOG"
SyncScript="/home/lidar/SCRIPTS/SH/MH_WINDLIDAR_FIG_SYNC.SH"
SyncLog="/home/lidar/SCRIPTS/SH/MH_WINDLIDAR_FIG_SYNC.LOG"

# Scan ID definition
# add lists of all IDs that apply
ScanID={}
# any PPI scan:
ScanID['PPI']   = [20, 24, 34, 37, 38, 48, 49]
# any RHI scan:
ScanID['RHI']   = [26, 32]
# any line-of-sight measurements:
ScanID['LOS']   = [20, 30, 41, 42, 44, 45]
# 89.99 degrees elevation, 0 degrees azimuth:
ScanID['LOS90'] = [30]
# vertical staring, for PBL detection:
ScanID['PBL']   = [20]
# beta calibration scan (low elevation, homogeneous boundary layer):
ScanID['CAL']   = [34]
# VAD scan (full PPI, 0 to 360 degrees azimuth, for wind fit):
ScanID['VAD']   = [37, 48, 49]
# any low level scans (low elevation):
ScanID['LOW']   = [20, 34, 38, 49]


# attributes for netcdf creation
AttDict={"location"  : 'Mace Head Atmospheric Research Station',
        "latitude"   : 53.33,
        "longitude"  : -9.9,
        "altitude"   : 20,
        "institute"  : 'NUI Galway, Ireland',
        "instrument" : 'WindCube 200S, Leosphere'
        }


# switches
SWITCH_REMOVE_BG = False    # remove background from plot (True), or plot background (False)
SWITCH_ZOOM      = False    # zoom in to background noise (change color bar limits, only for CNR) (True), or uses limits given in VarDict (False)
SWITCH_NC        = True     # uses existing netcdf files if in data path (True), or uses all text files in data path as input (False), or appends latest text file in data path to existing netcdf file in data path and removes this text file ('append')
SWITCH_OUTPUT    = True     # prints status messages on screen if run from command line (True)
SWITCH_TIMER     = True     # times the main processes while running the script, prints time elapsed since start of script if output is activated (True)
SWITCH_MODE      = 'LOW'    # calculates/plots only certain scan types ('VAD', 'LOW', 'LOS', 'LOS90'), or all scan types ('all')


# LOS zoom (range axis)
# plots the specified range only, if given; otherwise plots the whole profile
# dictionary keys are scan IDs, values are a list of plot limits in meter [min_range, max_range]
LOSzoom={}
LOSzoom[41] = [5500, 6500]
LOSzoom[42] = [2000, 3000]
LOSzoom[44] = [5000, 5500]
LOSzoom[45] = [2000, 6500]


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
                  "cols"  : ("time", "config_ID", "scan_ID", "LOS_ID", "azimuth", "elevation", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "names" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "m/s", "m/s", "dB", "no unit", "percent", "no unit"),
                  "longs" : ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "radial wind speed", "standard deviation of the radial wind speed", "carrier-to-noise ratio", "confidence index", "mean error", "status"),
                  "lims"  : ([], [], [], [], [-180,180], [0,90], [0,8000], [-20,20], [-1,1], [-30,10], [0,100], [], [0,1]),
                  "ty"    : ('d', 'i', 'i', 'i', 'd', 'd', 'i', 'd', 'd', 'd', 'd', 'd', 'b')
            },
        "VAD"  : {"fend"  : "_whole_radial",
                  "N"     : 2,
                  "cols"  : ("time", "range", 'speed', 'vertical', 'direction', 'confidence_index', 'rsquared', 'number_of_function_calls'),
                  "names" : ("time", "range", 'wspeed', 'w', 'wdir', 'confidence_index', 'rsquared', 'number_of_function_calls'),
                  "units" : ("seconds since 1970-01-01 00:00:00", "meter above ground level", "m s-1", "m s-1", "degree", "no unit", "no unit", "no unit"),
                  "longs" : ("time", "distance from sensor to center of each range gates along the line of sight", "wind_speed", "upward_air_velocity", "wind_from_direction", "confidence index", "R^2 of fit", "number of function calls for least spuare fit"),
                  "lims"  : ([], [0,8000], [0,25], [-2,2], [0,360], [0,100], [0,1], []),
                  "ty"    : ('d', 'i', 'd', 'd', 'd', 'd', 'd', 'b')
            },
        "cnr"  : {"fend"  : "_whole_radial",
                  "N"     : 9,
                  "cols"  : ("time", "config_ID", "scan_ID", "LOS_ID", "azimuth", "elevation", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "names" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "radial_wind_speed", "dev_radial_wind_speed", "CNR", "confidence_index", "mean_error", "status"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "m/s", "m/s", "dB", "no unit", "percent", "no unit"),
                  "longs" : ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "radial wind speed", "standard deviation of the radial wind speed", "carrier-to-noise ratio", "confidence index", "mean error", "status"),
                  "lims"  : ([], [], [], [], [-180,180], [0,90], [0,8000], [-20,20], [-1,1], [-30,10], [0,100], [], [0,1]),
                  "ty"    : ('d', 'i', 'i', 'i', 'd', 'd', 'i', 'd', 'd', 'd', 'd', 'd', 'b')
            },
        "beta" : {"fend"  : "_radial_beta",
                  "N"     : 7,
                  "cols"  : ("time", "config_ID", "scan_ID", "LOS_ID", "azimuth", "elevation", "range", "att_rel_beta"),
                  "names" : ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "att_rel_beta"),
                  "units" : ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "1/(m*sr)"),
                  "longs" : ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "attenuated relative backscatter"),
                  "lims"  : ([], [], [], [], [-180,180], [0,90], [0,8000], [-8,-3]),
                  "ty"    : ('d', 'i', 'i', 'i', 'd', 'd', 'i', 'd')
            },
        "spectra":{"fend" : "_spectra",
                   "N"    : 7,
                   "cols" : ("time", "config_ID", "scan_ID", "LOS_ID", "azimuth", "elevation", "range", "spectra"),
                   "names": ("time", "config_ID", "scan_ID", "LOS_ID", "azi", "ele", "range", "spectra"),
                   "units": ("seconds since 1970-01-01 00:00:00", "no unit", "no unit", "no unit", "degree", "degree", "m", "a.u."),
                   "longs": ("time", "configuration ID", "scan ID", "line-of-sight ID", "sensor_azimuth_angle", "sensor_elevation_angle", "distance from sensor to center of each range gates along the line of sight", "spectral density"),
                   "lims" : ([], [], [], [], [-180,180], [0,90], [0,8000], [0,5000]),
                   "ty"   : ('d', 'i', 'i', 'i', 'd', 'd', 'i', 'd')
            }
        }
