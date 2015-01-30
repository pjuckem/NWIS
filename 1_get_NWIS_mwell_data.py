__author__ = 'pfjuckem'
'''
retrieve a list of NWIS wells by county.  Keep wells with more than X msmts.

INPUT:
model domain    --> shapefile of model domain (polygon)
counties        --> shapefile of counties (polygon)
min_msmts       --> minimum number of measurements for keeping the well for futher analysis

OUTPUT:
state_lookup_full -->  file with site information for wells in NWIS that meet criteria


'''
import os
import arcpy
import geopandas as gpd
import numpy as np
import NWIS_utilities as nwispy

# Inputs
model_domain = 'D:/PFJData2/Projects/NAQWA/Cycle3/FWP/ARC/PEST/FluxTargetHUCS.shp'
#UTM83Z16_ft = {u'proj':u'utm', u'zone':16, u'datum':u'NAD83', u'units':u'us-ft', u'no_defs':True} #UTM83 zone 16 feet, manually defined
counties = 'D:/ARC/Basemaps/National/Boundaries/UScounties.shp'
#wi_counties_crs = 102973
#wi_counties_crs = crs.from_epsg(wi_counties_crs)
#inputEPSG = 4326 # lat long

startdate = '1970-01-01'
enddate = '2012-12-31'
#stateFIPS = '055'

#pump_gdf = gpd.GeoDataFrame(pumpdf, crs=UTM83Z16_ft)

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

if __name__ == '__main__':
    fipslist = open(countylistfile, 'w')
    arcpy.env.overwriteOutput = True
    arcpy.Intersect_analysis([model_domain, counties], fwpcounties, 'ALL')
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
    '''
    for state in stateFIPS:
        for cCounty in countyFIPS:
            sfront, sback = os.path.splitext(siteSCfile)
            sitefile = sfront + '_' + str(state) + str(cCounty) + sback
            dtwfront, dtwback = os.path.splitext(dtwSCfile)
            dtwfile = dtwfront + '_' + str(state) + str(cCounty) + dtwback
            nwispy.puller_by_state_county('site',state,cCounty,startdate,enddate,sitefile,'county')
            nwispy.puller_by_state_county('gwlevels',state,cCounty,startdate,enddate,dtwfile,'county')
            #nwispy.puller_by_latlong('site',stateFIPS,cCounty,startdate,enddate,sitefile,'county')
            #nwispy.puller_by_latlong('gwlevels',stateFIPS,cCounty,startdate,enddate,dtwfile,'county')
            statecounty = str(state+cCounty)
            fipslist.write(statecounty + ', ')
'''




'''
# read in info from daily values sites file
header_text = open(NWIS_daily_values_sites_file).readlines()
columns, header_rows = ft.NWIS_header(header_text)

df = pd.read_csv(NWIS_daily_values_sites_file, sep='\t', names=columns, skiprows=header_rows)

for n in df.site_no:

    # first test if the site is in the model domain
    site_location = Point(df[df.site_no == n][['dec_long_va', 'dec_lat_va']].get_values()[0])

    if site_location.within(bounds):
        text = ft.get_nwis(n, '00060')
        ofp = open(os.path.join(output_folder, '{}.txt'.format(n)), 'w')
        [ofp.write(line) for line in text]
        ofp.close()
'''