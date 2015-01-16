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
import geopandas as gpd
import pandas as pd

# Input files
infofile='D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/mwell_siteinfo_137.dat'
levelsfile='D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/mwell_dtw_137.dat'

# Outfiles
pdffile='extended_records.pdf'

# Settings
mode='GFLOW' # GFLOW or MODFLOW; writes either a tp file, or .hob file for MF2k observation process

min_msmts = 2
# A set of dictionaries to evaluate the accuracy of individual msmts.  Will need to interpret lev_errdict in light of the others.
# see: http://pubs.usgs.gov/of/2004/1238/pdf/gwcoding_Sect2-3.pdf
lev_errdict = {0:1, 1:0.1, 2:0.01, 9:10} # direct estimate of accuracy
src_errdict = {'A':0.01,'D':1,'G':0.01,'L':0.1,'M':10,'O':1,'R':1,'S':0.01,'Z':1} # agency msmts likely accurate; Drillers, owners, other likely from time of drilling
meth_errdict = {'A':1,'B':0.01,'C':0.01,'E':5,'F':0.01,'G':0.01,'H':0.01,'L':0.1,'M':0.01,'N':0.1,'R':1,'S':01,'T':0.01,'U':5,'V':0.01,'Z':1} # method accuracy.
# list of level_status_cds to discard, e.g. if well was being pumped (using ALL codes)
discard_lev_status_cd=['A','B','D','E','F','G','H','I','J','M','N','O','P','R','S','T','V','W','X','Z']
nmsmt_weightdict = {1:20, 10:5, 11:1}  # 1 msmt: var=20ft, 2-10 msmts: var=5ft, >=11msmts: var=1 ft

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

info_header=getheader(infofile,'agency_cd','\t')
levels_header=getheader(levelsfile,'agency_cd','\t')

info=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None)[1:]
levelsdata=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]

# convert to dataframes and set site_no as the index (as a long integer)
siteinfodf = pd.DataFrame(info)
siteinfodf.site_no = siteinfodf.site_no.astype(long)
siteinfodf.index = siteinfodf['site_no']
levelsdf = pd.DataFrame(levelsdata)
levelsdf['site_no'] = levelsdata['site_no'].astype(long)
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

5. Bin according to:
    A: water level accuracy; assign numeric value (1, 0.1, 0.01 ft)
    C: altitude accuracy; assign numeric value based on coded value (10, 1, 0.01 ft)
    C: source of measurement; assign numeric value based on an estimated relationship (fabricated)
    D: measurement method; assign numeric value based on an estimated relationship (fabricated)
    ...Sum up A - D = estimate of reliability.  Look at some output and evaluate how to convert this to a weight for pest
4. For remaining wells, add lat-long, then convert to UTM-ft;  compute average WL and std (& plot?)
6. process well ID, location, well top& bot, ave WL, and weight for writing out to MOD2OBS.
7. any additional processing for mod2obs.
'''
# Filter the dataset, evaluate accuracy codes, compute average WLs
# 1. if alt_acy_va is 0.01, keep it as it was surveyed by USGS
# 2. if msmt has a lev_status code, discard it.
# 3. count number of msmts per well.  If alt_acy_va <= 0.01 OR n_msmts > min msmts, then keep it.
# 4. Use the worst accuracy (using dictionaries above) of the lev_src_cd, lev_meth_cd, and lev_acy_cd fields.
# 5. Add accuracy from #4 to alt_acy_va
# 6. Limit to msmts within a specified time period
# 7. Compute average of all WL msmts per well and add as a column to siteinfodf
# 8. Compute weight for this target value as a function of line #5 and the number of msmts, with more msmts = better, as per dict above
levelsdf = levelsdf.copy().join(siteinfodf.loc[:, ['alt_acy_va']])  # join the desired column based on the indexes (siteID)
levelsdf = levelsdf.loc[~levelsdf['lev_status_cd'].isin(discard_lev_status_cd), :] # keep msmts w/out errors. Note: the "~" inside of the first "[" means "is NOT in"
levelsdf['n_msmts'] = levelsdf.groupby(['site_no'])['lev_va'].count()  # number wl msmts per well
# edit next line to keep wells w/ >minmsmt OR alt_accuracy = 0.01
levelsdf = levelsdf.loc[levelsdf['n_msmts'] >=min_msmts, :] # keep wells w/ more than the minimum # msmts

# if source isn't agency or
#levelsdf['lev_err'] =
levelsdf['cum_err'] = levelsdf.groupby(['site_no'])['lev_err'].sum()  # number wl msmts per well






'''
code not used, but useful for future reference:

siteobjcols = list(siteinfodf.select_dtypes(include=[object])) # list of colums w/ dtype = object (str or mix)
levelsobjcols = list(levelsdf.select_dtypes(include=[object])) # list of colums w/ dtype = object (str or mix)

'''