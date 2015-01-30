'''
program to QC groundwater levels from NWIS database and generate mod2obs file.

Required Inputs:
list of counties with NWIS Site Files and level files generated by "1_get_NWIS_mwell_data.py
required NWIS fields: site_no, alt_va, alt_acy_va, lev_dt, lev_va, lev_status_cd, well_depth_va, ....others

Output:  Mod2Obs files + PEST INS file and portion of PST file with head obs.
'''

import numpy as np
import datetime as dt
import geopandas as gpd
import pandas as pd
from fiona import crs
from shapely.geometry import Point, LineString, Polygon
import sys
import os
import flopy
import re

if sys.platform == 'darwin' or 'nix' in sys.platform:
    newline = '\r\n'
else:
    newline = '\n'

'''
SETUP SECTION
'''
# Outfiles
outpath = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/'
outshape = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC/PEST/targets/NWIS_targets.shp'

# Input files
calib_area = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC/PEST/FluxTargetHUCS.shp'
mfpath = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/working'
mfnam = 'FWPvert_uniformall.nam'
MFdis = 'FWPvert.dis'
origin = '826111.706, 15791723.564'  # origin or lower left of model grid in world coordinates of the model (UTM-ft)
uprleft = '916863.802, 16717285.046'
datadir = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/'
infofileprefix = 'mwell_siteinfo_'
levelfileprefix = 'mwell_dtw_'
countylistfile = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/countylist.dat'
#countylist = ['137']  # note: 141 = wood county = some wells in;some out

# MOD2OBS input:  *.spc and *.fig exported from GWvistas.  Need to update for multi-layer models.
driverfile = 'mod2obs_FWP.in'
bore_coords_file = 'FWP_bore_coords.crd'
bore_list_file = 'FWP_bore_list.crd'
well_sample_file = 'FWP_GWSI.smp'
hds_file = 'FWPvert_uniformall.hds'
insfile = 'FWP_GWSI.ins'
pest_head_obs_sect_file = 'FWP_pst_GWSI.dat'
spc_file = 'FWPvert.spc'
flow_or_trans = 'f'
inactive_thresh = 1.0e+30
time_units = 'd'
nlay = 1
extrapolation_limit = 1.0e+10
output_file = 'FWP_uniformall_GWSI.smp'
ss_arbitrary_date = '12/20/2012'
ss_arbitrary_time = '12:00:00'

# Settings
glac_only = True
start = dt.datetime(1970, 1, 1)  # start and end dates for evaluating WL data
end = dt.datetime(2013, 12, 31)
min_msmts = 2  # minimum number of wl msmts in NWIS to be included as a target (eliminates pre-1980 WCRs)
onlysubft = True  # another filter to retain only reasonably decent wells.
# A set of dictionaries to evaluate the accuracy of individual msmts. values are in feet.
# see: http://pubs.usgs.gov/of/2004/1238/
lev_errdict = {0:1.01, 1:0.1, 2:0.01, 9:10}  # direct estimate of msmt accuracy.
src_errdict = {'A':0.01,'D':5,'G':0.01,'L':0.1,'M':10,'O':5,'R':.1,'S':0.01,'Z':.1}  # agency msmts likely accurate; Drillers, owners, other likely from time of drilling
meth_errdict = {'A':1.01,'B':0.01,'C':0.01,'E':5,'F':0.01,'G':0.01,'H':0.01,'L':0.1,'M':0.01,'N':0.1,'R':1.01,'S':0.01,'T':0.01,'U':.1,'V':0.01,'Z':1.01}  # method accuracy.
# list of level_status_cds to discard. EG: if well was recently pumped, discard the msmt (but not well).(using ALL codes)
discard_lev_status_cd=['A','B','D','E','F','G','H','I','J','M','N','O','P','R','S','T','V','W','X','Z']
# list of coordinate_codes to keep. Keeping only wells with locations better or equal to 5 degree-seconds (F)
#keep_coord_cd = {'H', '1', '5', 'S', 'R', 'F'}  # might want to add 'T' (10 degree-seconds, or about 900ft latitude)
# keep_coords_cd now removed b/c many long-term wells are coded as 1-minute accuracy, but is likely better than that.
# house keeping
driverfile = os.path.join(outpath, driverfile)
bore_coords_file = os.path.join(outpath, bore_coords_file)
bore_list_file = os.path.join(outpath, bore_list_file)
well_sample_file = os.path.join(outpath, well_sample_file)
hds_file = os.path.join(outpath, hds_file)
insfile = os.path.join(outpath, insfile)
pest_head_obs_sect_file = os.path.join(outpath, pest_head_obs_sect_file)
spc_file = os.path.join(outpath, spc_file)
output_file = os.path.join(outpath, output_file)

'''
FUNCTIONS
'''
# Functions from: http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2':
            >> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def angle(orig, point1, point2, unit='deg'):
    """  Input is 2 or 3 x,y coordinate pairs supplied as arrays. """
    v1 = point1 - orig
    v2 = point2 - orig
    rad_angle = np.array([angle_between(v1, v2)])
    if unit == 'deg':
        deg_angle = np.degrees(rad_angle)
        return deg_angle
    elif unit == 'rad':
        return rad_angle
    else:
        print 'Error.  Please specify "deg" for degrees or "rad" for radians'

def xy2np(string):
    """
    Convert a string of x and y into a numpy array vector of x & y
    """
    x,y = string.split(',')
    x = float(x)
    y = float(y)
    xy = np.array([x, y])
    return xy

def calc_angle(str_orig, string1, string2=None, unit='deg'):
    """
    Wrap-up all of the above into a one-liner.
    If string2 is not provided, assume point2 is on y-axis and desired angle is from pos y-axis (from north).
    """
    orig = xy2np(str_orig)
    point1 = xy2np(string1)
    if string2:
        point2 = xy2np(string2)
    else:
        point2 = orig + np.array([0,10])
    arr_angle = angle(orig, point1, point2, unit)
    return arr_angle # return array

# get header info
def getheader(filename,indicator,delimiter):
    headervar=0
    infile=open(filename,'r').readlines()
    for line in infile:
        cline=line.strip().split(delimiter)
        if cline[0]==indicator:
            break
        else:
            headervar+=1
    return(headervar)


'''
MAIN CODE
'''

print "getting well info, water levels, and coordinates..."
countyfile = open(countylistfile, 'r')
countylist = countyfile.readline().split(',')
countylist = ''.join(countylist).split()  # convert to str, remove whitespaces
n = 0
for c in countylist:
    infofile = os.path.join(datadir, (infofileprefix + c + '.dat'))
    levelsfile = os.path.join(datadir, (levelfileprefix + c + '.dat'))
    
    info_header=getheader(infofile,'agency_cd','\t')
    levels_header=getheader(levelsfile,'agency_cd','\t')
    
    if n==0:
        info=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None, comments='garbalygook')[1:]  # don't use # as default comment -- it's in some well names
        siteinfodf = pd.DataFrame(info)
        levelsdata=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]
        levelsdf = pd.DataFrame(levelsdata)
    else:
        info1=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None, comments='garbalygook')[1:]
        siteinfo1 = pd.DataFrame(info1)
        siteinfodf = siteinfodf.append(siteinfo1)
        levelsdata1=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]
        levels1 = pd.DataFrame(levelsdata1)
        levelsdf = levelsdf.append(levels1)
    n += 1

# convert to dataframes and set site_no as the index
siteinfodf.site_no = siteinfodf.site_no.astype(long)
siteinfodf.index = siteinfodf['site_no']
levelsdf['site_no'] = levelsdf['site_no'].astype(long)
levelsdf.index = levelsdf['site_no']

# Clean-up fields and change dtype so can aggrigate
#  This function doesn't work; not clear why. Needs more work. Commands work when not part of a function...
def clean_to_floats(df, fieldlist):
    for field in fieldlist:
        df = df.loc[df[field] != ('\D+'), :] # remove rows with non-float values.
        df = df.loc[df[field] != (''), :] # remove empty rows (ideas to add this to line above?).
        df[field] = df[field].astype(float)
#siteselectcols = ['dec_lat_va', 'dec_long_va', 'alt_va', 'alt_acy_va', 'well_depth_va']
#levelsselectcols = ['lev_va', 'lev_acy_cd']
#clean_to_floats(siteinfodf, siteselectcols)
#clean_to_floats(levelsdf, levelsselectcols)

siteinfodf = siteinfodf.loc[siteinfodf['well_depth_va']!=('\D+'), :]  # remove rows with non-float values.
siteinfodf = siteinfodf.loc[siteinfodf['well_depth_va']!=(''), :]  # don't keep blank rows
siteinfodf.well_depth_va = siteinfodf.well_depth_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['dec_lat_va']!=('\D+'), :]
siteinfodf = siteinfodf.loc[siteinfodf['dec_lat_va']!=(''), :]
siteinfodf.dec_lat_va = siteinfodf.dec_lat_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['dec_long_va']!=('\D+'), :]
siteinfodf = siteinfodf.loc[siteinfodf['dec_long_va']!=(''), :]
siteinfodf.dec_long_va = siteinfodf.dec_long_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['alt_va']!=('\D+'), :]
siteinfodf = siteinfodf.loc[siteinfodf['alt_va']!=(''), :]
siteinfodf.alt_va = siteinfodf.alt_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['alt_acy_va']!=('\D+'), :]
siteinfodf = siteinfodf.loc[siteinfodf['alt_acy_va']!=(''), :]
siteinfodf.alt_acy_va = siteinfodf.alt_acy_va.astype(float)
levelsdf = levelsdf.loc[levelsdf['lev_va']!=('\D+'), :]
levelsdf = levelsdf.loc[levelsdf['lev_va']!=(''), :]
levelsdf.lev_va = levelsdf.lev_va.astype(float)
levelsdf = levelsdf.loc[levelsdf['lev_acy_cd']!=('\D+'), :]
levelsdf = levelsdf.loc[levelsdf['lev_acy_cd']!=(''), :]
levelsdf.lev_acy_cd = levelsdf.lev_acy_cd.astype(float)

# Filter the dataset, evaluate accuracy codes, compute average WLs
levelsdf = levelsdf.copy().join(siteinfodf.loc[:, ['alt_acy_va']])  # join based on the indexes
levelsdf = levelsdf.loc[~levelsdf['lev_status_cd'].isin(discard_lev_status_cd), :] # keep msmts w/out errors. Note: the "~"
# add columns then decode the dictionaries
levelsdf['src_err'] = levelsdf['lev_src_cd'].copy()
levelsdf['meth_err'] = levelsdf['lev_meth_cd'].copy()
levelsdf['levacy_err'] = levelsdf['lev_acy_cd'].copy()
levelsdf = levelsdf.replace({'src_err':src_errdict}).copy()
levelsdf = levelsdf.replace({'meth_err':meth_errdict}).copy()
levelsdf = levelsdf.replace({'levacy_err':lev_errdict}).copy()
levelsdf['cum_err'] = levelsdf['levacy_err'] + levelsdf['meth_err'] + levelsdf['src_err']

if onlysubft:
    print "limiting wells to those with cumulative measurement errors < 1 ft"
    levelsdf = levelsdf[(levelsdf['cum_err'] < 1.0)]  # This retains many decent USGS wells.
levelsdf['n_msmts'] = levelsdf.groupby(['site_no'])['lev_va'].count()
levelsdf = levelsdf[(levelsdf['alt_acy_va'] <= 0.01) | (levelsdf['n_msmts'] >= min_msmts)]  # min threshold for keeping

# Limit by dates, process for errors, compute ave WL.
levelsdf['lev_dt'] = pd.to_datetime(levelsdf['lev_dt'])
levelsdf = levelsdf[(levelsdf['lev_dt'] >= start) & (levelsdf['lev_dt'] <= end)]
levelsdf['first'] = levelsdf.groupby(['site_no'])['lev_dt'].min()
levelsdf['last'] = levelsdf.groupby(['site_no'])['lev_dt'].max()
tdiff = levelsdf['last'] - levelsdf['first']  # minimum of 1 day
tdiff = tdiff.reset_index()
tdiff.columns = ['site_no', 'diff']
tdiff.index = tdiff['site_no']
#tdiff['td'] = [1.157407 * 10**-14 * float(td) for td in tdiff['diff']]  # convert to float; convert to days.
# This worked with NP v1.7 & PD v0.14, but broke after update to NP1.8 & PD0.15. wouldn't work when tried to revert.
# Temporary hack while Pandas team works to better handle timedelta dtypes....
tdiff['td'] = [float(str(td).split()[0]) for td in tdiff['diff']]  # Note: units are now presumed to be days rather than nanoseconds.
tdiff2 = tdiff.groupby(tdiff.index).agg(np.mean)  # warning, the column 'site_no' gets changed
tdiff2 = tdiff2.drop('site_no', axis=1)  # line above altered 'site_no' column (not index), so remove it for safety.

lev2 = levelsdf.groupby('site_no').agg(np.mean)  # compute mean WL (all fields), then aggregates by site number
lev2 = lev2.copy().join(tdiff2.loc[:])  # add time differences

# Assign Weights based on altitude err, cumulative msmt err, # msmts, and time span of msmts
# Duration of msmts & number of msmts deemed most important for weights.  Secondarily, highly accurate altitudes are
# likely associated with accurate locations (RTK-GPS), so alt_err also a priority.
# Priority for weighting: 1) time span, 2) # measurements, 3) alt accuracy, 4) cum msmt err
print "assigning weights"
lev2['weight'] = -1.0
lev2['group'] = 'NoGroup'
for i, s in lev2.iterrows():
    #   map or GPS location?    Measured with tape or estimated?   approx 2 msmts per year?    > 20 yrs of record?
    if (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  #BEST of the best
        lev2.loc[i, 'weight'] = 50
        lev2.loc[i, 'group'] = 'long-t_mwell'
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  # better of the best
        lev2.loc[i, 'weight'] = 40
        lev2.loc[i, 'group'] = 'long-t_mwell'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  # best
        lev2.loc[i, 'weight'] = 29
        lev2.loc[i, 'group'] = 'long-t_mwell'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  # best
        lev2.loc[i, 'weight'] = 24
        lev2.loc[i, 'group'] = 'long-t_mwell'

    #   map or GPS location?    Measured with tape or estimated?  approx 1 msmt per year?    5 - 19.9 yrs of record?
    elif (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # best
        lev2.loc[i, 'weight'] = 30
        lev2.loc[i, 'group'] = 'mid-t_mwell'
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # best
        lev2.loc[i, 'weight'] = 25
        lev2.loc[i, 'group'] = 'mid-t_mwell'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # good
        lev2.loc[i, 'weight'] = 17
        lev2.loc[i, 'group'] = 'mid-t_mwell'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # good
        lev2.loc[i, 'weight'] = 16
        lev2.loc[i, 'group'] = 'mid-t_mwell'

    #   map or GPS location?    Measured with tape or estimated?   approx 1 msmt per season?    > 1 yr of record?
    elif (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # good
        lev2.loc[i, 'weight'] = 20
        lev2.loc[i, 'group'] = 'shrt-t_well'
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # good
        lev2.loc[i, 'weight'] = 15
        lev2.loc[i, 'group'] = 'shrt-t_well'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # fair
        lev2.loc[i, 'weight'] = 11
        lev2.loc[i, 'group'] = 'shrt-t_well'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # fair
        lev2.loc[i, 'weight'] = 10
        lev2.loc[i, 'group'] = 'shrt-t_well'

    # Note: earlier filters should have removed wells with <min_msmts unless they had really good locational accuracy
    #   map or GPS location?    Measured with tape or estimated?   take what we can get..........
    elif (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.td[i] < 365.242):  # fair
        lev2.loc[i, 'weight'] = 9
        lev2.loc[i, 'group'] = 'few-msmt_well'
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.td[i] < 365.242):  # poor
        lev2.loc[i, 'weight'] = 8
        lev2.loc[i, 'group'] = 'few-msmt_well'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.td[i] < 365.242):  # poor
        lev2.loc[i, 'weight'] = 5
        lev2.loc[i, 'group'] = 'few-msmt_well'
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.td[i] < 365.242):  # poor
        lev2.loc[i, 'weight'] = 4
        lev2.loc[i, 'group'] = 'few-msmt_well'

    else: lev2.loc[i, 'weight'] = -2.0  # Should ID problems

# Join-in site information
lev3 = lev2.join(siteinfodf.loc[:, ['station_nm', 'county_cd', 'dec_lat_va', 'dec_long_va', 'coord_acy_cd',
                                     'nat_aqfr_cd', 'aqfr_cd', 'aqfr_type_cd', 'alt_va', 'well_depth_va']])

# Keep only glacial wells w/ accuracy of 5 seconds or better (BR can tolerate less accuracy due to flow sys detail)
if glac_only:
    print "removing non-glacial wells"
    lev3 = lev3[((lev3['aqfr_cd'] == '100SDGV') | (lev3['nat_aqfr_cd'] == 'N100GLCIAL') | (lev3['aqfr_cd'] == '') |
                  (lev3['aqfr_cd'] == '110QRNR') | (lev3['aqfr_cd'] == '110SDGVP') | (lev3['aqfr_cd'] == '110SANDP'))]
    #lev3 = lev3.loc[lev3['coord_acy_cd'].isin(keep_coord_cd), :]  # keep sites with location w/in +/- X seconds (F=5sec)

print "re-projecting coordinates to UTM-ft, and clipping to the calibration area"
x, y = lev3['dec_long_va'], lev3['dec_lat_va']
xy = zip(x, y)
points = gpd.GeoSeries([Point(x, y) for x, y in xy])
lev3 = lev3.reset_index()  # have to revert to a basic index to get points to match up
lev3['geometry'] = points
inputEPSG = 4326 # lat long
crs = crs.from_epsg(inputEPSG)
lev3 = gpd.GeoDataFrame(lev3, crs=crs)
UTM83Z16_ft = {u'proj':u'utm', u'zone':16, u'datum':u'NAD83', u'units':u'us-ft', u'no_defs':True}  # UTM83 zone 16 feet
lev3.to_crs(crs=UTM83Z16_ft, inplace=True)

# Use Geopandas to determine which points are in/out.
# Hack to process the returned bool series (there has to be a better way...)
# Note: can't use arcpy.clip because object dtype columns are dropped, such as station name.
calib = gpd.read_file(calib_area)
geom = calib.geometry
keep = [geom.contains(pt) for pt in lev3['geometry']]
lst = []
for r, ser in enumerate(keep):
    val = str(keep[r]).split()[1]
    if val == 'True': val = True
    elif val == 'False': val = False
    else:
        print "error"
        break
    lst.append(val)
lev3.loc[:,'keep'] = lst
lev3 = lev3.loc[lev3['keep']==True, :]  # remove wells outside the calibration area from the dataframe
lev3 = lev3.set_index('site_no')  # reset index after converting coordinates and clipping

print "Truncating names to 10-spaces for PEST"
labels = []
names = lev3['station_nm'].tolist()
for name in names:
    last = name.split('/')[-1]
    try:
        last[10]  # if longer than 10 characters, truncate
        last = last[-10:]
    except:
        pass
    last = last.strip()
    labels.append(last)
lev3.loc[:,'labels'] = labels

print "     Checking for duplicate names and fixing them"
lev3['dup'] = lev3.duplicated(subset='labels')
dups = lev3.copy()[lev3['dup'] == True]
n = dups['dup'].count()
while n > 0:
    newlabels = []
    names = dups['labels'].tolist()
    for i, name in enumerate(names):
        last = name[-7:]  # take the last 7 characters
        last = str(i+1) + '_' + last  # add an index
        newlabels.append(last)
    dups.loc[:,'labels'] = newlabels
    lev3.loc[:,'newlabels'] = dups.loc[:,'labels']
    lev3['labels'] = np.where(lev3.loc[:, 'newlabels'].isnull(), lev3['labels'], lev3['newlabels'])  # .loc for bool series
    lev3['dup'] = lev3.duplicated(subset='labels')
    dups = lev3.copy()[lev3['dup'] == True]
    n = dups['dup'].count()

print "assigning targets to layers"
m = flopy.modflow.Modflow(model_ws=mfpath)
nf = flopy.utils.mfreadnam.parsenamefile(os.path.join(mfpath, mfnam), {})
dis = flopy.modflow.ModflowDis.load(os.path.join(mfpath, MFdis), m, nf)
tops = np.zeros((dis.nrow, dis.ncol))
tops[:, :] = dis.top.array
bots = np.zeros((dis.nlay, dis.nrow, dis.ncol))
bots[:, :, :] = dis.botm.array
nrow = dis.nrow
ncol = dis.ncol

lev3['wellelev'] = lev3['alt_va'] - lev3['well_depth_va']
wellbotelev = np.array(lev3['wellelev'].tolist())
rad_angle = calc_angle(origin, uprleft, unit='rad')  # radian_angle from true north
geom = lev3['geometry'].tolist()
x_orig = float(origin.split(',')[0])
y_orig = float(origin.split(',')[1])
xcoords = [geom[i].x for i,j in enumerate(geom)]
ycoords = [geom[i].y for i,j in enumerate(geom)]
#  compute row,col indexes for a rotated grid
wellcolsidx = [int(np.floor((x - x_orig - (np.tan(rad_angle) * (ycoords[i] - y_orig)))/(1000.0 / np.cos(rad_angle)))) for i, x in enumerate(xcoords)]  # columns are 1000ft
wellrowsidx = [nrow-1 - int(np.floor((y - y_orig + (np.tan(rad_angle) * (xcoords[i] - x_orig)))/(1000.0 / np.cos(rad_angle)))) for i, y in enumerate(ycoords)]

layer = np.full_like(wellbotelev, -999)  # -999 error flag
topatwells = np.array([tops[i,j] for i,j in zip(wellrowsidx, wellcolsidx)])
for k in range(bots.shape[0]):
    layKbotatwells = np.array([bots[k,i,j] for i,j in zip(wellrowsidx, wellcolsidx)])
    layer = np.where(wellbotelev >= layKbotatwells, k+1, layer)  # designate the lowest layer the target can be in
    layer = np.where(wellbotelev > topatwells, 999, layer)  # a flag if the well bot is above LS -- what's going on?
minval = layer.min()
if minval < 0:
    print "At least one target has a well bottom below the bottom of the model.  They will be removed as targets"
maxval = layer.max()
if maxval >= 999:
    print "At least one target has a well bottom above the land surface.  They have been assigned to layer 1," \
          "but they should be inspected using the shapefile (layer = 999) and site_information files."

lev3['layer'] = layer
lev3['x'] = xcoords
lev3['y'] = ycoords
lev3['date'] = ss_arbitrary_date
lev3['time'] = ss_arbitrary_time
lev3['target'] = lev3['alt_va'] - lev3['lev_va']
lev3 = lev3[(lev3['layer'] != -999)]  # removing wells below the model bottom

# Format for Mod2Obs, including unique 10-digit IDs
# Setup for mod2obs1; could use regular mod2obs if remove "a" values from ofp.write, below.
ofp = open(driverfile, 'w')
ofp.write('{1}{0}{2}{0}{3}{0}{4}{0}a{0}{5}{0}\
{6}{0}{7:<10.2e}{0}{8}{0}{9}{0}{10}{0}\
{11}{0}{12:<10.2e}{0}{13}{0}a{0}'.format(newline,
                                 spc_file,
                                 bore_coords_file,
                                 bore_list_file,
                                 well_sample_file,
                                 hds_file,
                                 flow_or_trans,
                                 inactive_thresh,
                                 time_units,
                                 ss_arbitrary_date,
                                 ss_arbitrary_time,
                                 nlay,
                                 extrapolation_limit,
                                 output_file))
ofp.close()

lev3.to_csv(bore_coords_file, columns=[u'labels', u'x', u'y', u'layer'], sep=' ', index='', header='')
lev3.to_csv(bore_list_file, columns=[u'labels'], sep=' ', index='', header='')
lev3.to_csv(well_sample_file, columns=[u'labels', u'date', u'time', u'target'],
            sep=' ', float_format='%.3f', index='', header='')
# make ins file
names_upper = [i.upper() for i in list(lev3['labels'])]
lev3['Name'] = names_upper
lev3['w'] = 'w'
lev3['hashname'] = '#' + lev3['Name'] + '#'
lev3['bangname'] = '!' + lev3['Name'] + '!'
## write the ins file
open(insfile, 'w').write('pif #\n')
with open(insfile, 'a') as file:
    lev3.to_csv(file, columns=[u'hashname',u'w',u'w',u'w',u'bangname'],
                  sep=' ', index='', header='')
# make pst file parts
lev3.to_csv(pest_head_obs_sect_file, columns=[u'Name',u'target',u'weight',u'group'],
                  sep=' ', float_format='%.3f', index='', header='', line_terminator=newline)
print '\nSuccessful completion.  Mod2Obs files are located here:\n' \
      '{}'.format(outpath)
lev3 = lev3.drop(['keep', 'dup', 'newlabels', 'date', 'time', 'w', 'hashname', 'bangname'], axis=1)
lev3.to_file(outshape)