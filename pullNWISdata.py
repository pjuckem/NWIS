__author__ = 'pfjuckem'
'''
retrieve a list of NWIS wells by county.  Keep wells with more than X msmts.

INPUT:
calib_area      --> shapefile of calibration area (polygon)
counties        --> shapefile of counties (polygon)

OUTPUT:
state_lookup_full -->  file with site information for wells in NWIS that meet criteria
'''
import os
import arcpy
import geopandas as gpd
import numpy as np
import NWIS_utilities as nwispy

# Inputs
calib_area = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC/PEST/FluxTargetHUCS.shp'
counties = 'D:/ARC/Basemaps/National/Boundaries/UScounties.shp'

startdate = '1970-01-01'
enddate = '2012-12-31'

# outputs
output_folder = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup'
siteSCfile = 'mwell_siteinfo.dat'
siteSCfile = os.path.join(output_folder, siteSCfile)
dtwSCfile = 'mwell_dtw.dat'
dtwSCfile = os.path.join(output_folder, dtwSCfile)
arc_folder = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC'
fwpcounties = 'fwpcounties.shp'
fwpcounties = os.path.join(arc_folder, fwpcounties)
countylistfile = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/MODFLOW/PESTsetup/countylist.dat'

# ### MAIN ###
fipslist = open(countylistfile, 'w')
arcpy.env.overwriteOutput = True
arcpy.Intersect_analysis([calib_area, counties], fwpcounties, 'ALL')
domain_counties = gpd.read_file(fwpcounties)
allFIPS = domain_counties['FIPS'].values
# countyFIPS = domain_counties['CNTY_FIPS'].values
# stateFIPS = np.array(domain_counties['STATE_FIPS'].values)
# stateFIPS = np.unique(stateFIPS)

for fips in allFIPS:
    state = fips[:2]
    county = fips[-3:]
    sfront, sback = os.path.splitext(siteSCfile)
    sitefile = sfront + '_' + str(fips) + sback
    dtwfront, dtwback = os.path.splitext(dtwSCfile)
    dtwfile = dtwfront + '_' + str(fips) + dtwback
    nwispy.puller_by_state_county('site',state,county,startdate,enddate,sitefile,'county')
    nwispy.puller_by_state_county('gwlevels',state,county,startdate,enddate,dtwfile,'county')
    #nwispy.puller_by_latlong('site',stateFIPS,cCounty,startdate,enddate,sitefile,'county')
    #nwispy.puller_by_latlong('gwlevels',stateFIPS,cCounty,startdate,enddate,dtwfile,'county')
    fipslist.write(str(fips) + ', ')
fipslist.close()
