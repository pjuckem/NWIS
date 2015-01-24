# program to QC groundwater levels from NWIS database

# Required Inputs:
# NWIS Site File (with station info; tab delimited)
# NWIS groundwater levels file (with levels and dates; tab delimited)
# required NWIS fields: agency_cd, site_no, alt_va, alt_acy_va, lev_dt, lev_va, lev_status_cd, well_depth_va, qw_count_nu, reliability_cd
#
# Coordinates File (with coordinates in model coordinate system for each well)
# This can be created by exporting csv of stations and WGS84 coordinates into ArcMap
# A character needs to be added to ends of station numbers (e.g. '434238088592501n') so that Arc treats them as strings
# Coordinates exported from Arc should have columns site_no2,POINT_X,POINT_Y
#
# The sorting algorithm below might be cleaner and more transparent if reformulated in terms of cumulative error
# e.g. 5 ft. + 3.5 ft. for alt_accuracy of <=5 ft., and std deviation in measurements of 3.5
# would still have to come up with quantities to represent the lower perceived error for wells with Water Quality measurements, artesian, etc.

import numpy as np
from collections import defaultdict
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
import arcpy
import geopandas as gpd
import pandas as pd
from fiona import crs
from shapely.geometry import Point, LineString, Polygon
from scipy import stats
import os

# Input files
calib_area = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC/PEST/FluxTargetHUCS.shp'
datadir = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/'
infofileprefix = 'mwell_siteinfo_'
levelfileprefix = 'mwell_dtw_'
countylist = ['137']  # note: 141 = wood county = some wells in;some out
#infofile='D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/mwell_siteinfo_97.dat'
#levelsfile='D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/mwell_dtw_97.dat'

# Outfiles
pdffile='extended_records.pdf'
#levelsout = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/mwell_dtw_97.out'
outshape = 'D:/ARC/junk/NWIS_FWPtargets.shp'
outshapeclip = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC/PEST/targets/NWIS_FWPtargets.shp'

# Settings
arcpy.env.overwriteOutput = True
mode='GFLOW' # GFLOW or MODFLOW; writes either a tp file, or .hob file for MF2k observation process
glac_only = True
start = dt.datetime(1970, 1, 1)
end = dt.datetime(2013, 12, 31)
min_msmts = 2
onlysubft = True
# A set of dictionaries to evaluate the accuracy of individual msmts.  Will need to interpret lev_errdict in light of the others.
# see: http://pubs.usgs.gov/of/2004/1238/pdf/gwcoding_Sect2-3.pdf
lev_errdict = {0:1.01, 1:0.1, 2:0.01, 9:10} # direct estimate of accuracy
src_errdict = {'A':0.01,'D':5,'G':0.01,'L':0.1,'M':10,'O':5,'R':.1,'S':0.01,'Z':.1} # agency msmts likely accurate; Drillers, owners, other likely from time of drilling
meth_errdict = {'A':1.01,'B':0.01,'C':0.01,'E':5,'F':0.01,'G':0.01,'H':0.01,'L':0.1,'M':0.01,'N':0.1,'R':1.01,'S':0.01,'T':0.01,'U':.1,'V':0.01,'Z':1.01} # method accuracy.

# list of level_status_cds to discard, e.g. if well was being pumped (using ALL codes)
discard_lev_status_cd=['A','B','D','E','F','G','H','I','J','M','N','O','P','R','S','T','V','W','X','Z']
keep_coord_cd = {'H', '1', '5', 'S', 'R', 'F'} # might want to add 'T' (10 degree-seconds, or about 900ft latitude)
#nmsmt_weightdict = {1:20, 10:5, 11:1}  # 1 msmt: var=20ft, 2-10 msmts: var=5ft, >=11msmts: var=1 ft

print "getting well info, water levels, and coordinates..."

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

n = 0
for c in countylist:    
    infofile = os.path.join(datadir, infofileprefix + c + '.dat')
    levelsfile = os.path.join(datadir, levelfileprefix + c + '.dat')
    
    info_header=getheader(infofile,'agency_cd','\t')
    levels_header=getheader(levelsfile,'agency_cd','\t')
    
    if n==0:
        info=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None)[1:]
        siteinfodf = pd.DataFrame(info)
        levelsdata=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]
        levelsdf = pd.DataFrame(levelsdata)
    else:
        info1=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None)[1:]
        siteinfo1 = pd.DataFrame(info1)
        siteinfodf = siteinfodf.append(siteinfo1)
        levelsdata1=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]
        levels1 = pd.DataFrame(levelsdata1)
        levelsdf = levelsdf.append(levels1)
    n += 1
'''    

info_header=getheader(infofile,'agency_cd','\t')
levels_header=getheader(levelsfile,'agency_cd','\t')

info=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None)[1:]
levelsdata=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]
'''

# convert to dataframes and set site_no as the index (as a long integer)
#siteinfodf = pd.DataFrame(info)
siteinfodf.site_no = siteinfodf.site_no.astype(long)
siteinfodf.index = siteinfodf['site_no']
#levelsdf = pd.DataFrame(levelsdata)
levelsdf['site_no'] = levelsdf['site_no'].astype(long)
levelsdf.index = levelsdf['site_no']

#  This function doesn't work; not clear why. Commands work when not part of a function...
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
siteinfodf = siteinfodf.loc[siteinfodf['dec_lat_va']!=('\D+'), :]  # remove rows with non-float values.
siteinfodf = siteinfodf.loc[siteinfodf['dec_lat_va']!=(''), :]  # don't keep blank rows
siteinfodf.dec_lat_va = siteinfodf.dec_lat_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['dec_long_va']!=('\D+'), :]  # remove rows with non-float values.
siteinfodf = siteinfodf.loc[siteinfodf['dec_long_va']!=(''), :]  # don't keep blank rows
siteinfodf.dec_long_va = siteinfodf.dec_long_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['alt_va']!=('\D+'), :]  # remove rows with non-float values.
siteinfodf = siteinfodf.loc[siteinfodf['alt_va']!=(''), :]  # don't keep blank rows
siteinfodf.alt_va = siteinfodf.alt_va.astype(float)
siteinfodf = siteinfodf.loc[siteinfodf['alt_acy_va']!=('\D+'), :]  # remove rows with non-float values.
siteinfodf = siteinfodf.loc[siteinfodf['alt_acy_va']!=(''), :]  # don't keep blank rows
siteinfodf.alt_acy_va = siteinfodf.alt_acy_va.astype(float)
levelsdf = levelsdf.loc[levelsdf['lev_va']!=('\D+'), :]  # remove rows with non-float values.
levelsdf = levelsdf.loc[levelsdf['lev_va']!=(''), :]  # don't keep blank rows
levelsdf.lev_va = levelsdf.lev_va.astype(float)
levelsdf = levelsdf.loc[levelsdf['lev_acy_cd']!=('\D+'), :]  # remove rows with non-float values.
levelsdf = levelsdf.loc[levelsdf['lev_acy_cd']!=(''), :]  # don't keep blank rows
levelsdf.lev_acy_cd = levelsdf.lev_acy_cd.astype(float)


'''
4. For remaining wells, add lat-long, then convert to UTM-ft;  compute average WL and std (& plot?)
6. process well ID, location, well top& bot, ave WL, and weight for writing out to MOD2OBS.
7. any additional processing for mod2obs.
'''
# Filter the dataset, evaluate accuracy codes, compute average WLs

levelsdf = levelsdf.copy().join(siteinfodf.loc[:, ['alt_acy_va']])  # join the desired column based on the indexes (siteID)
levelsdf = levelsdf.loc[~levelsdf['lev_status_cd'].isin(discard_lev_status_cd), :] # keep msmts w/out errors. Note: the "~" inside of the first "[" means "is NOT in"
levelsdf['n_msmts'] = levelsdf.groupby(['site_no'])['lev_va'].count()  # number wl msmts per well
levelsdf = levelsdf[(levelsdf['alt_acy_va'] <= 0.01) | (levelsdf['n_msmts'] >= min_msmts)]

# add columns then decode the dictionaries
levelsdf['src_err'] = levelsdf['lev_src_cd'].copy()
levelsdf['meth_err'] = levelsdf['lev_meth_cd'].copy()
levelsdf['levacy_err'] = levelsdf['lev_acy_cd'].copy()
levelsdf = levelsdf.replace({'src_err':src_errdict}).copy()
levelsdf = levelsdf.replace({'meth_err':meth_errdict}).copy()
levelsdf = levelsdf.replace({'levacy_err':lev_errdict}).copy()
#levelsdf['worst_err'] = levelsdf[['levacy_err', 'meth_err', 'src_err']].max(axis=1) # msmt can be no more accurate than the worst of these 3 accuracy indicators
levelsdf['cum_err'] = levelsdf['levacy_err'] + levelsdf['meth_err'] + levelsdf['src_err']

# Limit by dates, process for errors, compute ave WL, add as column to site dataframe. 
levelsdf['lev_dt'] = pd.to_datetime(levelsdf['lev_dt'])
levelsdf = levelsdf[(levelsdf['lev_dt'] >= start) & (levelsdf['lev_dt'] <= end)]
levelsdf['first'] = levelsdf.groupby(['site_no'])['lev_dt'].min()
levelsdf['last'] = levelsdf.groupby(['site_no'])['lev_dt'].max()
tdiff = levelsdf['last'] - levelsdf['first']  # minimum of 1 day
tdiff = tdiff.reset_index()
tdiff.columns = ['site_no', 'diff']
tdiff.index = tdiff['site_no']
#tdiff['td'] = [1.157407 * 10**-14 * float(td) for td in tdiff['diff']]  # convert to float; convert to days.  NOTE: this worked with NP v1.7 & PD v0.14, but broke after update to NP1.8 & PD0.15....& wouldn't work when tried to revert.  WTF?
# Temporary hack while Pandas team works to better handle timedelta dtypes....
tdiff['td'] = [float(str(td).split()[0]) for td in tdiff['diff']]  # Note: units are now presumed to be days rather than nanoseconds.
tdiff2 = tdiff.groupby(tdiff.index).agg(np.mean)  # warning, the column 'site_no' gets changed
tdiff2 = tdiff2.drop('site_no', axis=1)  # line above altered 'site_no' column (not index), so remove it.

# Had planned to ID the max errors and eliminate them if they're not the primary source of data.  NONE OF THESE ARE WORKING!!  Would need to MultiIndex or use panels, which is for the future....
#  Instead:  What if simply remove all msmts where cum_err >= 1.  Thus, keep only the best msmts.  Turns out that most msmts are good, as long as consider meth_acy_cd = U (unknown) a good msmt (often ignored recently)
if onlysubft:
    #levelsdf = levelsdf[(levelsdf['worst_err'] < 1.0)]  # keep only msmts with worst error less than 1 foot. This retains many decent USGS wells.
    levelsdf = levelsdf[(levelsdf['cum_err'] < 1.0)]  # keep only msmts with worst error less than 1 foot. This retains many decent USGS wells.
# Not sure I like this approach of computing and removing outliers -- not ruthless enough.  Plus, isn't grouping by both siteID and worst_err.  Why not simply remove the worst msmts instead?.  Then only keep the best
#levelsdf['mean_err'] = levelsdf.groupby(['site_no'])['worst_err'].mean()  # Review Note: This is not implemented correctly
#levelsdf['std_err'] = levelsdf.groupby(['site_no'])['worst_err'].std()  # Review Note: This is not implemented correctly
#levelsdf = levelsdf[np.abs(levelsdf.worst_err - levelsdf.mean_err) <= (3 * levelsdf.std_err)]  # keep only rows for which the "worst_error" is not an outlier

levelsdf['n_msmts'] = levelsdf.groupby(['site_no'])['lev_va'].count()  # Re-count number wl msmts per well before use in computations

lev2 = levelsdf.groupby('site_no').agg(np.mean)  # computes mean WL (& all fields), then aggregates by site number
lev2 = lev2.copy().join(tdiff2.loc[:])  # add time differences

# Assign Weights based on altitude err, cumulative msmt err, # msmts, and time span of msmts
# Given challenges with pin-point accuracy (DEM, MF cell size), # msmts and time span seem more important for weights
# Secondarily, highly accurate altitudes are likely associated with accurate locations (RTK-GPS, so alt_err also a priority
# Priority for weighting: 1) time span, 2) # measurements, 3) alt accuracy, 4) cum msmt err
lev2['weight'] = -1.0
for i, s in lev2.iterrows():
    #   map or GPS location?    Measured with tape or estimated?   approx 2 msmts per year?    > 20 yrs of record?
    if (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  #BEST of the best
        lev2.loc[i, 'weight'] = 50
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  # better of the best
        lev2.loc[i, 'weight'] = 40
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  # best
        lev2.loc[i, 'weight'] = 21
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] >= 40 and lev2.td[i] >= (20* 365.242)):  # best
        lev2.loc[i, 'weight'] = 20

    #   map or GPS location?    Measured with tape or estimated?  approx 1 msmt per year?    5 - 19.9 yrs of record?
    elif (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # best
        lev2.loc[i, 'weight'] = 30
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # best
        lev2.loc[i, 'weight'] = 25
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # good
        lev2.loc[i, 'weight'] = 17
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] >= 5 and lev2.td[i] >= (5 * 365.242)):  # good
        lev2.loc[i, 'weight'] = 16

    #   map or GPS location?    Measured with tape or estimated?   approx 1 msmt per season?    > 1 yr of record?
    elif (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # good
        lev2.loc[i, 'weight'] = 20
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # good
        lev2.loc[i, 'weight'] = 15
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # fair
        lev2.loc[i, 'weight'] = 11
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] >= 4 and lev2.td[i] >= 365.242):  # fair
        lev2.loc[i, 'weight'] = 10

    #   map or GPS location?    Measured with tape or estimated?   take what we can get..........
    # Note: earlier filters should have removed wells with <min_msmts unless they had really good locational accuracy
    elif (lev2.alt_acy_va[i] <= .5 and lev2.cum_err[i] <= .1 and lev2.n_msmts[i] < 4 and lev2.td[i] < 365.242):  # fair
        lev2.loc[i, 'weight'] = 9
    elif (lev2.alt_acy_va[i] <= 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] < 4 and lev2.td[i] < 365.242):  # poor
        lev2.loc[i, 'weight'] = 8
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] <= 1.0 and lev2.n_msmts[i] < 4 and lev2.td[i] < 365.242):  # poor
        lev2.loc[i, 'weight'] = 5
    elif (lev2.alt_acy_va[i] > 1.0 and lev2.cum_err[i] > 1.0 and lev2.n_msmts[i] < 4 and lev2.td[i] < 365.242):  # poor
        lev2.loc[i, 'weight'] = 4

    else: lev2.loc[i, 'weight'] = -2.0  # Should ID problems

# Join-in site information
lev3 = lev2.join(siteinfodf.loc[:, ['station_nm', 'county_cd', 'dec_lat_va', 'dec_long_va', 'coord_acy_cd',
                                     'nat_aqfr_cd', 'aqfr_cd', 'aqfr_type_cd', 'well_depth_va']])
'''
Need to deal with the fact that some well_depth_va values are actually the elevation of the well bottom!!
Will simply need to check:  If welldepthva > 1000, don't compute well elevation as DEM-welldepthva (just leave it)

Or:  If site_tp_cd = GW-TH (test hole, not completed as a well), don't use it.
'''

# Keep only glacial wells w/ accuracy of 5 seconds or better (BR can tolerate less accuracy due to flow sys detail)
if glac_only:
    lev3 = lev3[(lev3['aqfr_cd'] == '100SDGV') | (lev3['nat_aqfr_cd'] == 'N100GLCIAL')]
    lev3 = lev3.loc[lev3['coord_acy_cd'].isin(keep_coord_cd), :]  # keep sites with location w/in +/- X seconds (F=5seconds)

# generate a 'geometry' field in order to create a geospatial dataframe and shapefile
x, y = lev3['dec_long_va'], lev3['dec_lat_va']
xy = zip(x, y)
points = gpd.GeoSeries([Point(x, y) for x, y in xy])
lev3 = lev3.reset_index()  # have to revert to a basic index to get points to match up

lev3['geometry'] = points
inputEPSG = 4326 # lat long
crs = crs.from_epsg(inputEPSG)
lev3 = gpd.GeoDataFrame(lev3, crs=crs)
UTM83Z16_ft = {u'proj':u'utm', u'zone':16, u'datum':u'NAD83', u'units':u'us-ft', u'no_defs':True} #UTM83 zone 16 feet, manually defined
#pump_gdf = gpd.GeoDataFrame(pumpdf, crs=UTM83Z16_ft)
#lev3.to_crs(epsg=26915, inplace=True) # UTM83 zone 15
lev3.to_crs(crs=UTM83Z16_ft, inplace=True) # UTM83 zone 16, feet

# Use Geopandas to determine which points are in/out.  Hack to process the returned series (there has to be a better way...)
# Note: couldn't use arcpy.clip because object dtype columns were dropped, such as station name.
calib = gpd.read_file(calib_area)
geom = calib.geometry
keep = [geom.contains(pt) for pt in lev3['geometry']]
ref = lev3.index.tolist() # get the indexes
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

print "Truncating names to 10-spaces for PEST"
labels = []
names = lev3['station_nm'].tolist()
for name in names:
    last = name.split('/')[-1]
    try:
        last[10]  # is this string longer than 10 characters?  If so, truncate
        last = last[-10:]
    except:
        pass
    labels.append(last)
lev3.loc[:,'labels'] = labels
print "Checking for duplicate names and fixing"
lev3 = lev3.set_index('site_no')
lev3['dup'] = lev3.duplicated(subset='labels')
dups = lev3.copy()[lev3['dup'] == True]
n = dups['dup'].count()
if n > 0:
    newlabels = []
    names = dups['labels'].tolist()
    #indxs = dups.index()
    for i, name in enumerate(names):
        last = name.split('-')[-1]
        last = last[-7:]  # regardless, take the last 7 characters
        last = str(i+1) + '-' + last  # add an index
        newlabels.append(last)
    dups.loc[:,'labels'] = newlabels


    #lev3.loc[:,'labels'] = dups.loc[:,'labels']
    # try an np.where statement instead
    # lev3.where(dups.loc[:], dups.loc['labels'], lev3.loc['labels']
'''
# Fix duplicate names
dupscount=0
print "modifying any duplicate names by choosing new 10-digit strings from site numbers..."
print "(duplicate names occur because only a 10-digit subset of the 15-digit ID is used for names;\nsome programs like GFLOW or PEST have observation name length limits)"

for name in names.itervalues():
    nameslist.append(name)
num_dups=len(nameslist)-len(np.unique(np.array(nameslist)))
i=0
while num_dups>0:
    dupscount=0
    i+=1
    wells2delete=[]
    for well in names.iterkeys():
        if len(names[well])==0:
            wells2delete.append(well)
            continue
        count=nameslist.count(names[well])
        if count>1:
            dupscount+=1
            oldname=names[well]
            if i<=4:
                names[well]=well[4-i:(-1-i)]+oldname[-5:]
                print names[well]
            else:
                names[well]=well[i:]+oldname[-5:]
    print "fixed %s duplicate names!" %(dupscount)
    for well in wells2delete:
        print 'deleting well: %s\n' %(well)
        del names[well]
    nameslist=[]
    for name in names.itervalues():
        nameslist.append(name)
    num_dups=len(nameslist)-len(np.unique(np.array(nameslist)))
    if num_dups>0:
        print "still duplicates; trying new 10-digit strings from site numbers.."
    else:
        break
print "done fixing duplicates!"
'''
'''
# format for Mod2Obs
if mode=='mod2obs':
    ofp=open('NWIS_FWP_mod2obs.csv','w')
    print 'NWIS_FWP_mod2obs.csv'
    ofp.write('Name,POINT_X,POINT_Y,WL,ScreenTop,ScreenBot\n')
    for category in ['best','good','fair','poor']:
        for well in wells:
            if category in names[well]:
                if WellDepth_elev[well]==None:
                    continue
                else:
                    ofp.write('%s,%s,%s,%s,%s,%s\n' %(names[well],coords[well][0],coords[well][1],levels[well],WellDepth_elev[well]+screen_length,WellDepth_elev[well]))
    ofp.close()
print "finished OK"
'''
# Format for Mod2Obs, including unique 10-digit IDs
# 1. Use misc code to convert lat long to UTMft
# 2. Clip to area of interest
# 3. Format targets for Mod2Obs

#levelsdf.to_csv(levelsout)
'''
code not used, but useful for future reference:

siteobjcols = list(siteinfodf.select_dtypes(include=[object])) # list of colums w/ dtype = object (str or mix)
levelsobjcols = list(levelsdf.select_dtypes(include=[object])) # list of colums w/ dtype = object (str or mix)


for reference:
multis = [i for i in y if y[i] > 1]

aq_wells['MNW2_top'] = np.where(aq_wells[top_open] == 1, (disobj.top[(aq_wells.row - 1), (aq_wells.column - 1)]), aq_wells[top_open])
df.loc[(df['BBB'] > 25) | (df['CCC'] >= 75), 'AAA'] = 0.1
df.ix[df.AAA >= 5,'BBB'] = -1
df.replace({"col1": di})
'''