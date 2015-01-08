# set of functions to use web services to get data from NWISweb
# to parse the resulting text (rdb) files, and to plot results
#
# original by Mike Fienen, additions by Howard Reeves
#

import urllib
import re
import os
import numpy as np
from datetime import datetime
from datetime import date
from datetime import timedelta
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 # this strange bit of code makes fonts 
                                  # editable in Adobe Illustator





mpl.rcParams['font.sans-serif']          = 'Univers 57 Condensed' #'Arial'
mpl.rcParams['font.serif']               = 'Times'
mpl.rcParams['font.cursive']             = 'Zapf Chancery'
mpl.rcParams['font.fantasy']             = 'Comic Sans MS'
mpl.rcParams['font.monospace']           = 'Courier New'
mpl.rcParams['pdf.compression']          = 0
mpl.rcParams['pdf.fonttype']             = 42

ticksize = 6
mpl.rcParams['legend.fontsize']  = 6
mpl.rcParams['axes.labelsize']   = 8
mpl.rcParams['xtick.labelsize']  = ticksize
mpl.rcParams['ytick.labelsize']  = ticksize

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from collections import defaultdict
import math

# ############
# function definitions
# ############

#--Function to get data from NWISweb using web services with list of station numbers
# original from Mike Fienen, 'datatype' added to allow for pulling other web services
#           datatype ------  meaning
#             dv             daily values, pulled in original by Mike
#            site            sitefile information, with output set to expanded,
#                            gets file that includes aquifer code, altitude of well, etc.
#           gwlevels         rdb file of manual groundwater level measurements (tapedowns)
# HWR, August 15, 2012
def puller_by_stations(datatype,cStations,sttime,endtime,outfilename):
    print 'Pulling GWSI data for list of stations'
    urlParts = ['http://waterservices.usgs.gov/nwis/TYPE/?format=rdb,1.0&sites=',
                '&startDT=',
                '&endDT=',
                '&siteType=GW']
    
    fullURL = re.sub('TYPE',datatype,urlParts[0])
    for i,cs in enumerate(cStations):
        fullURL += cs
        if i+1 < len(cStations):
            fullURL += ','
    if re.match('site',datatype):
        fullURL+='&siteOutput=expanded'
    fullURL += urlParts[1] + sttime + urlParts[2] + endtime + urlParts[3]
    print 'Using URL:\n%s' %(fullURL)
    datastream = urllib.urlopen(fullURL).read()
    open(outfilename,'wb').write(datastream)
    print 'File download complete'
    print 'output written to %s' %(outfilename)
#--Function to get data from NWISweb using web services with state and county specified
# original from Mike Fienen, 'datatype' added to allow for pulling other web services
#           datatype ------  meaning
#             dv             daily values, pulled in original by Mike
#            site            sitefile information, with output set to expanded,
#                            gets file that includes aquifer code, altitude of well, etc.
#           gwlevels         rdb file of manual groundwater level measurements (tapedowns)
# HWR, August 15, 2012
def puller_by_state_county(datatype,cState,cCounty,sttime,endtime,outfilename,stat_county):
    if stat_county == 'county':
        print 'Pulling data for %s in %s' %(cCounty, cState)
    elif stat_county == 'state':
        print 'Pulling data for %s' %(cState)
    state_lookup_full,state_lookup_abbrev,county_lookup,state_full_abbrev_lookup = read_state_county_FIPS()    
    stateFIPS,countyFIPS = get_county_and_state_FIPS(cState,cCounty,state_lookup_full,state_lookup_abbrev,county_lookup)
    urlParts = ['http://waterservices.usgs.gov/nwis/TYPE/?format=rdb',
            '&countyCd=',
            '&stateCd=',
            '&startDT=',
            '&endDT=',
            '&siteType=GW',]
    fullURL = re.sub('TYPE',datatype,urlParts[0])
    if re.match('site',datatype):
        fullURL+='&siteOutput=expanded'
    if stat_county == 'county':
        fullURL += urlParts[1] + stateFIPS + countyFIPS + urlParts[3] + sttime + urlParts[4] + endtime + urlParts[5]
    elif stat_county == 'state':
        if len(cState) > 2:
            cState = state_lookup_abbrev[cState]
        fullURL += urlParts[2] + cState + urlParts[3] + sttime + urlParts[4] + endtime + urlParts[5]
        
    print 'Using URL:\n%s' %(fullURL)
    datastream = urllib.urlopen(fullURL).read()
    open(outfilename,'wb').write(datastream)
    print 'File download complete'    
    print 'output written to %s' %(outfilename)

#--Function to make a lookup of State and County FIPS codes
def read_state_county_FIPS():
    '''
    read_state_county_FIPS()
    A function to read FIPS files for lookup and return two kinds of output.
    Mike Fienen - 7/16/2012
    <mnfienen *at* usgs *dot* gov>
    
    INPUT:
    none --> the variables are all already named
    
    OUTPUT:
    state_lookup_full -->  dictionary with keys of full state names 
                    and elements of state codes
    state_lookup_abbrev --> dictionary with keys of state name 
                    abbreviations and elements of state codes
    county_lookup --> a 3 column numpy array of strings with county codes,
                    state codes, and county names
    '''
    # filenames for the State and County lookups
    state_file = os.path.join('..','NWIS_meta_data','STATE_FIPS.csv')
    county_file = os.path.join('..','NWIS_meta_data','COUNTY_FIPS.csv')
    #
    # read in the county file --> READ THE "HARD WAY"
    #
    # make lists to hold the data we will read
    state_codes = []
    county_codes = []
    county_names = []
    
    indat = open(county_file,'r').readlines()  # read the whole file into memory as a list
    headers = indat.pop(0)                     # remove the 0th line and put into headers
    for line in indat:                         # go through the remaining file line by line
        tmp = line.strip().split(',')          # remove newline chars and split on ','
        state_codes.append(tmp[1])             # leave as strings to keep zero padding
        county_codes.append(tmp[2])
        county_names.append(tmp[3].lower())    # leave the names as strings but !make lower case!
    # convert the lists to numpy arrays--usefule later to enable using np.nonzero
    state_codes = np.atleast_1d(state_codes)
    county_codes = np.atleast_1d(county_codes)
    county_names = np.atleast_1d(county_names)
        
    # now mash them up into one numpy array
    county_lookup = np.hstack((state_codes,county_codes,county_names))
    # reshape it to be columns again and transpose the results (.T)
    county_lookup = county_lookup.reshape(3,len(county_lookup)/3).T
    
    #
    # read in the state file --> USES NP.GENFROMTXT
    #
    state_vals = np.genfromtxt(state_file,dtype=None,names=True,delimiter=',')
    # force the state names and abbreviations to lower case for later comparisons
    fullnames = state_vals['State_Name']
    for i in np.arange(len(fullnames)):
        fullnames[i] = fullnames[i].lower()
    abbrevnames = state_vals['State_Abbreviation']
    for i in np.arange(len(abbrevnames)):
        abbrevnames[i] = abbrevnames[i].lower()
    # now make the output lookup dictionaries to return
    state_lookup_full = dict(zip(fullnames,
                                 state_vals['FIPS_Code']))
    state_lookup_abbrev = dict(zip(abbrevnames,
                                   state_vals['FIPS_Code']))
    state_abbrev_full = dict(zip(fullnames,abbrevnames))
    return state_lookup_full,state_lookup_abbrev,county_lookup,state_abbrev_full

#--Function to lookup state and county FIPS codes and associate them with state and county
def get_county_and_state_FIPS(cState,cCounty,state_lookup_full,
                              state_lookup_abbrev,county_lookup):
    '''
    get_county_and_state_FIPS(cState,cCounty,state_lookup_full,
                              state_lookup_abbrev,county_lookup)
    A function to lookup a state and county set of FIPS codes using 
    string representations of the states and counties, or optionally
    a FIPS for state and string for county. 
    N.B.--> some of the counties are actually, parishes, etc., so some further
                        QC may be in order!
    Mike Fienen - 7/16/2012
    <mnfienen *at* usgs *dot* gov>
    
    INPUT:
    cState --> state to lookup, can be FIPS code, full name, or 2 letter abbrev.
    cCounty --> County name (due to note above, must include descriptor such as "county" etc.
    state_lookup_full -->  dictionary with keys of full state names 
                    and elements of state codes
    state_lookup_abbrev --> dictionary with keys of state name 
                    abbreviations and elements of state codes
    county_lookup --> a 3 column numpy array of strings with county codes,
                    state codes, and county names
    OUTPUT:
    cFIPS --> Dictionary with keys of "state" and "county" and elements with the FIPS codes
    '''
    #
    # first be flexible about how to handle the state code ambiguity
    #
    try:  # is it in integer?
        state_FIPS = int(cState)
    except: # if not, see, by length, if it's an abbreviation or full name
        cState = cState.lower()  # force to lower case
        lenstate = len(cState)
        if lenstate == 2:
            try:
                state_FIPS = state_lookup_abbrev[cState]
            except KeyError:
                raise KeyError('State not found')
        else:
            try:
                state_FIPS = state_lookup_full[cState]
            except KeyError:
                raise KeyError( 'State not found')
    state_FIPS = str(state_FIPS).zfill(2) # convert back to a string padded with zeros on the left
    #
    # next, find a match in county_lookup for both state FIPS code and county name
    # returning the county FIPS code
    #
    try: # is it already an integer?
        county_FIPS = int(cCounty)
    except: # if not, use the string to look up
        cCounty = cCounty.lower()  # first force to lower case for the comparison        
        # now, find the indices within county_lookup that match both county name and state FIPS code
        indies = (np.where(county_lookup[:,0]==state_FIPS) and np.where(county_lookup[:,2]==cCounty))[0]
        if len(indies)>1:
            print "ambiguous county and state match"
        else:
            county_FIPS = county_lookup[indies,1][0]
    county_FIPS=str(county_FIPS).zfill(3) # convert back to string and pad with zeros on left
    return state_FIPS,county_FIPS
        
        
#--Function to read a file from NWIS query output
def NWIS_dv_reader(infile):
    '''
    NWIS_dv_reader(infile)
    A function to read in an NWIS file generated using USGS webservices.
    Mike Fienen - 7/16/2012
    <mnfienen *at* usgs *dot* gov>

    renamed by HWR to _dv_reader to distinguish it from readers for
    site file and tapedown data  8/15/2012
    
    INPUT:
    infile --> the name of an input file in USGS RDB (tab-delimited) format
    
    OUTPUT:
    indat --> a dictionary with keys corresponding to site numbers and each
            element being a dictionary with keys date and depth to water.
    '''
    # tell the user what's happening
    print 'Reading daily values from file: %s' %(infile)
    #
    # set up a couple initial variables
    #
    # format for reading the date --> formats noted at 
    #        http://docs.python.org/library/datetime.html (bottom of the page)
    indatefmt = "%Y-%m-%d"
    
    
    # open the text file and read all lines into a variable "tmpdat"
    tmpdat = open(infile,'r').readlines()
    
    # make empty lists to temporarily hold the data
    Site_ID = []   # site ID
    dates = []     # date of measurement (daily values)
    DTW = []       # depth to water below land surface (feet)
    prov_code = [] # provisional code: [P] is provisional, [A] is accepted
    
    # loop over the input data, keep only proper data rows. Parse and assign to lists
    for lnum, line in enumerate(tmpdat):
        # first read the lookup information from the header of the file
        if ("data for the following" in line.lower()):
            nWells = int(re.findall("[0-9]+",line)[0])
            statnums = []
            countynums = []
            for cwell in np.arange(nWells):
                nextline = lnum+1+cwell
                tmp = tmpdat[nextline].strip().split()
                statnums.append(tmp[2])
                countynums.append(tmp[3])
            station_lookup = dict(zip(statnums,countynums))
                
        if (('usgs' in line.lower()) and ('#' not in line)):
            tmp = line.strip().split() # strip newline off the end and split on whitespace
            Site_ID.append(tmp[1])
            dates.append(datetime.strptime(tmp[2],indatefmt)) #convert date to a time tuple
            DTW.append(tmp[3])
            prov_code.append(tmp[4].lower()) # --> note conversion to lower case!


    Site_ID = np.array(Site_ID,dtype=str) # convert to a numpy array
                                          # keep as str to not lose trailing zeroes (!)
    DTW = np.array(DTW,dtype=float) # convert to numpy array as a float containing the measurements
    dates = np.array(dates,dtype=object) #convert to numpy array as object since these are time tuples
    prov_code = np.array(prov_code,dtype=str) # convert to a numpy array as string 
    unique_sites = np.unique(Site_ID) # get a unique list of the sites
    
    # make an empty dictionary to include the results
    indat = dict()
    # now loop over the unique sites and parse the results by site number
    for csite in unique_sites:
        cindices = np.nonzero(Site_ID == csite)[0] # find the indices for the current site
        indat[csite] = {'dates':dates[cindices],
                        'DTW':DTW[cindices],
                        'prov_code':prov_code[cindices]}
    # tell the user we are done!
    print "Reading %s complete!" %(infile)
    return indat,station_lookup

#--Function to read the sitefile information read by puller_  functions
#  with datatype='site'
def NWIS_sitefile_reader(infile):
    '''
    NWIS_sitefile_reader(infile)

    A function to read the txt file of site information generated
    by the 'puller' function with datatype='site'
    Reads latitude, longitude, altitude of land surface at well,
    national aquifer code, local aquifer code, site name, well depth.
    Returns a dictionary with the stationid as the key and
    these values as a secondary dictionay (defaultdict(dict)).

    Howard W. Reeves - 8/15/2012
    <hwreeves *at* usgs *dot* gov>

    INPUT:
    infile --> filename generated by puller with datatype='site'

    OUTPUT:
    siteinfo --> defaultdict(dict),
        primary key=stationid
        secondary key=GWSI code
            site name,station_nm      
            latitude, dec_lat_va      
            longitude, dec_long_va
            altitude, alt_va
            datum, alt_datum_cd
            local aquifer code, aqfr_cd
            national aquifer code, nat_aqfr_cd
            well depth, well_depth_va
    '''
    
    siteinfo=defaultdict(dict)
    targetcodes={'site_no':1,
                 'station_nm':1,
                 'dec_lat_va':1,
                 'dec_long_va':1,
                 'alt_va':1,
                 'alt_datum_cd':1,
                 'aqfr_cd':1,
                 'nat_aqfr_cd':1,
                 'well_depth_va':1,
                 }
    targetcolumn=dict()

    print 'Reading sitefile information from file: %s' %(infile)

    #open the file and read until line doesn't start with #
    IN=open(infile,'r')
    header=1
    while header==1:
        line=IN.readline()
        if not re.match('^#',line.strip()):
            header=0
    #find column in tab-delimited file with desired header value
    vals=re.split('\t',line)
    column=0
    for item in vals:
        if item in targetcodes:
            targetcolumn[item]=column
        column+=1
    #read last header line
    line=IN.readline()
    #now read tab-delimited information and put into siteinfo
    for line in IN:
        vals=re.split('\t',line.strip())
        stationid=vals[targetcolumn['site_no']]
        for item in targetcodes:
            if not re.match('site_no',item):
                siteinfo[stationid][item]=vals[targetcolumn[item]]
    IN.close()
    return siteinfo

#--Function to read the gwlevels read by puller_  functions
#  with datatype='gwlevels'
def NWIS_gwlevels_reader(infile):
    '''
    NWIS_gwlevels_reader(infile)

    A function to read the txt file of groundwater levels generated
    by the 'puller' function with datatype='gwlevels'
    Returns dictionary with key=station number and numpy arrays of date
    and values 

    Howard W. Reeves - 8/15/2012
    <hwreeves *at* usgs *dot* gov>

    INPUT:
    infile --> filename generated by puller with datatype='gwlevels'

    OUTPUT:
    gwlevels --> dictionary with values read from file
    '''
    gwlevels=dict()
    targetcodes={'site_no':1,
                 'lev_dt':1,
                 'lev_tm':1,
                 'lev_va':1,
                 }
    targetcolumn=dict()

    print 'Reading gwlevels from file: %s' %(infile)

    #open the file and read until line doesn't start with #
    INDAT=open(infile,'r')
    header=1
    while header==1:
        line=INDAT.readline()
        if not re.match('^#',line.strip()):
            header=0
    #find column in tab-delimited file with desired header value
    vals=re.split('\t',line)
    column=0
    for item in vals:
        if item in targetcodes:
            targetcolumn[item]=column
        column+=1
    #read last header line
    line=INDAT.readline()
    #now read tab-delimited information and put into gwlevels
    gw=defaultdict(list)
    dates=defaultdict(list)
    for line in INDAT:
        vals=re.split('\t',line.strip())
        stationid=vals[targetcolumn['site_no']]
        if re.match('^\d',vals[targetcolumn['lev_va']].strip()):
            day=vals[targetcolumn['lev_dt']]
            tim=vals[targetcolumn['lev_tm']]
            if re.search('0',day):
                (year,month,dt)=re.split('-',day.strip())
                if re.match('^\d',tim.strip()):
                    (hr,mn)=re.split(':',tim.strip())
                else:
                    (hr,mn)=('12','0')
                when=datetime(int(year),int(month),int(dt),int(hr),int(mn))
                dates[stationid].append(when)
                gw[stationid].append(vals[targetcolumn['lev_va']])


    for station in gw.iterkeys():
        npgw = np.array(gw[station],dtype=float) # convert to numpy array as a float containing the measurements
        npdates = np.array(dates[station],dtype=object) #convert to numpy array as object since these are time tuples
        gwlevels[station]={'dates':npdates,
                           'levels':npgw,}
        
    return gwlevels
        
#--Function to plot time series in segments if there is a gap of over a given
#  number of days so that gaps are not filled with a straight line
def time_series_segments(dates,DTW,indices,symbol,gapsize):
    '''
    time_series_segments(time,DTW,indices,symbol,gapsize)

    A function to plot a time series in segments if there is a gap of over a given
    number of days so that gaps are not filled with a straight line. Figure
    and axes already set up in matplotlib in the calling routine NWIS_plotter.
    Set as a function because there are two time series: approved and provisional

    Howard W. Reeves - 8/16/2012
    <hwreeves *at* usgs *dot* gov>

    INPUT:
    dates --> vector of datetime values
    DTW --> depth to water value to be plotted
    indices --> indices to be plotted, vector set in NWIS_plotter
    symbol --> desired symbol for plot
    gapsize --> length in days to split segment

    OUTPUT:
    sends information to axes defined in NWIS_plotter
    '''
    nseg=1
    endindex=[]
    for i in range(1,len(dates[indices])):
        deltat=dates[indices][i]-dates[indices][i-1]
        if deltat.days > gapsize:
            nseg+=1
            endindex.append(i)
    if nseg==1:
        plt.plot(dates[indices],DTW[indices],symbol)
    else:
        srt=0
        endindex.append(len(dates[indices]))
        for i in range(0,nseg):
            plt.plot(dates[indices[srt:endindex[i]]],DTW[indices[srt:endindex[i]]],symbol)
            srt=endindex[i]
    return

#--Function to compute moving average of NWIS time series
def NWIS_moving_average(dates,DTW,tapedwn,datestp):
    '''
    NWIS_moving_average(dates,DTW,tapedwn,datestp)

    A function to compute a 5-day moving average for depth to water.
    Also computes standard deviation of DTW for each day of the water
    year and returns vectors for plotting.

    Howard W. Reeves - 8/16/2012
    <hwreeves *at* usgs *dot* gov>

    INPUT:
    dates --> vector of datetime values
    DTW --> vector of depth to water values
    tapedwn --> vector of tapedown depth to water levels
    datestp --> vector of dates for tapedowns

    OUTPUT:
    daynp --> numpy array of day of water year
    avenp --> numpy array of average DTW for each day of water year
    avemnp --> average minus standard deviation
    avepnp --> average plus standard deviation
    currday --> array of datetime days of water year for current water year
    currentnp --> numpy array of DTW for each day of current water year
    tapeday --> array of datetime days of water year for tapedowns in current year
    tapecurrent --> numpy array of tapedown DTW for current water year
    datesvec --> array of datetime days of a water year

    note -  dates in water year set into 2000 just for plotting, in the
            plotting function only months are marked.  This allows use
            of matplotlib date plotting features rather than having to
            build the ticks based on the day of water year vector.
    '''
    #make a default dict that breaks dates into water year and day of water year
    wateryearmin=3000
    wateryearmax=0
    depthtowater=defaultdict(dict)
    wydate=defaultdict(dict)
    thisyear=date.today().year
    thismonth=date.today().month
    if thismonth < 10:
        thiswateryear=thisyear
    else:
        thiswateryear=thisyear+1
    currentyes=0
    currday=[]
    currentnp=[]
    oneday=timedelta(days=1)
    for i in range(0,len(dates)):
        daterecord=dates[i]
        wateryear=daterecord.year
        month=daterecord.month
        if month >= 10:
            wateryear=wateryear+1
        if wateryear < wateryearmin:
            wateryearmin=wateryear
        if wateryear > wateryearmax:
            wateryearmax=wateryear
        firstdayyear=wateryear-1
        firstdayday=1
        firstdaymonth=10
        firstdate=datetime(firstdayyear,firstdaymonth,firstdayday)
        delta=daterecord-firstdate
        dayofwateryear=delta.days
        depthtowater[dayofwateryear][wateryear]=DTW[i]
    dayvec=[]
    avevec=[]
    avepvec=[]
    avemvec=[]
    currentvec=[]
    currentx=[]
    datesvec=[]
    firstday=datetime(2000,10,1)
    for day in range(3,364):
        summ=0
        ssumm=0
        pnts=0
        for i in range(wateryearmin,wateryearmax+1):
            for j in range(-2,3):
                inquestion=day+j
                if inquestion in depthtowater:
                    if i in depthtowater[inquestion]:
                        summ=summ+depthtowater[inquestion][i]
                        ssumm=ssumm+depthtowater[inquestion][i]**2
                        pnts+=1
        if pnts==0:
            continue
        else:
            datesvec.append(firstday+oneday*day)
            average=summ/pnts
            dayvec.append(day)
            avevec.append(average)
        if pnts<=1:
            pnts=1
            variance=(ssumm-(summ**2/pnts))/pnts
        else:
            variance=(ssumm-(summ**2/pnts))/(pnts-1)
        if variance <=0:
            stdev=0
        else:
            stdev=math.sqrt(variance)
        avepvec.append(average+stdev)
        avemvec.append(average-stdev)
        
        if day in depthtowater:
            if thiswateryear in depthtowater[day]:
                currentyes=1
                currentvec.append(depthtowater[day][thiswateryear])
                currday.append(firstday+oneday*day)
        #check to see if there are any tapedowns in the current water year
        tapeday=[]
        tapevec=[]
        tapenp=[]
        tapeyes=0
        if len(tapedwn)>0:
            for i in range(0,len(tapedwn)):
                if datestp[i].month<10:
                    wateryear=datestp[i].year
                else:
                    wateryear=datestp[i].year+1
                if wateryear == thiswateryear:
                    tempday=datestp[i].day
                    tempmonth=datestp[i].month
                    if month <10:
                        tempdate=datetime(2001,tempmonth,tempday)
                    else:
                        tempdate=datetime(2000,tempmonth,tempday)
                    tapeday.append(tempdate)
                    tapevec.append(tapedwn[i])
                    tapeyes=1
    #make and return numpy arrays
    
    daynp=np.array(dayvec,dtype=int)
    avenp=np.array(avevec,dtype=float)
    avemnp=np.array(avemvec,dtype=float)
    avepnp=np.array(avepvec,dtype=float)
    
    if currentyes==1:
        currentnp=np.array(currentvec,dtype=float)
    if tapeyes==1:
        tapenp=np.array(tapevec,dtype=float)
    
    return daynp,avenp,avemnp,avepnp,currentyes,currday,currentnp,datesvec,tapeyes,tapeday,tapenp

#--Function to plot moving average of depth to water from daily values
def NWIS_plot_moving_average(indat,gwdat,sitedat,plot_format,disp_plot=False):
    '''
    NWIS_plot_moving_average(indat,gwdat,sitedat,plot_format,disp_plot=False)

    Function to plot 5-day moving average of daily values downloaded
    from NWIS

    Howard W. Reeves - 8/16/2012
    <hwreeves *at* usgs *dot* gov>

    INPUT:
    indat --> a dictionary of data returned from puller 'dv'
    gwdat --> a dictionary of data returned from puller 'gwlevels'
    sitedat --> a dictionary of data returned from puller 'sites'
    plot_format --> format to save files of plots. can be '.png' or '.pdf'
    disp_plot --> a Boolean flag determining whether or not to display the plots onscreen
    
    OUTPUT:
    a set of files with plots of depth to water (reversing y axis) and
    named with station ID_movingave
    '''
    codeflag=0
    for cstation in sitedat.keys():
        #
        # first parse all the data for the current station
        #
        npcstation=np.array([cstation],dtype=str)
        if npcstation[0] in indat:
            DTW = indat[cstation]['DTW'] # get the depth to water for the current station
            datesdv = indat[cstation]['dates'] # pull the dates as datetime objects
            if cstation in gwlevels:
                tapedwn = gwdat[cstation]['levels']
                datestp = gwdat[cstation]['dates']
            else:
                tapedwn=[]
                datestp=[]
            # check if there are more than 1000 (3 years)
            # if so, take a chance on gaps and plot water year moving average
            if len(DTW)>1000:
                fig = plt.figure() 
                print 'computing moving average for {0:s}'.format(cstation)
                ax3=fig.add_subplot(211)
                ax3d=ax3.twinx()
                (daynp,avenp,avemnp,avepnp,currentyes,currdaynp,currentnp,datesvec,tapeyes,tapeday,tapenp)=NWIS_moving_average(datesdv,DTW,tapedwn,datestp)
                ymin=math.floor(avemnp.min()/5)*5
                ymax=math.ceil(avepnp.max()/5)*5
                ax3.plot(datesvec,avenp,'b-',datesvec,avepnp,'r-',datesvec,avemnp,'g-')
                # add line for current water year if values exist
                if currentyes==1:
                    if currentnp.min() < ymin:
                        ymin=math.floor(currentnp.min()/5)*5
                    if currentnp.max() > ymax:
                        ymax=math.ceil(currentnp.max()/5)*5
                    ax3.plot(currdaynp,currentnp,'k-')
                if tapeyes==1:
                    if tapenp.min() < ymin:
                        ymin=math.floor(tapenp.min()/5)*5
                    if tapenp.max() > ymax:
                        ymax=math.ceil(tapenp.max()/5)*5
                    ax3.plot(tapeday,tapenp,'bx')
                        
                ax3.set_ylim((ymax,ymin))
                altmin=float(sitedat[cstation]['alt_va'])-ymax
                altmax=float(sitedat[cstation]['alt_va'])-ymin
                ax3d.set_ylim([altmin,altmax])
                ax3d.set_ylabel('ALTITUDE OF WATER LEVEL IN FEET',fontsize=10)
                ax3.set_ylabel('DEPTH TO WATER IN FEET',fontsize=10)
                ax3.set_xlabel('WATER YEAR',fontsize=10)
                ax3.set_xlim((datetime(2000,10,1),datetime(2001,9,30)))
                monthloc=mdates.MonthLocator()
                monthfmt=mdates.DateFormatter('%b')
                ax3.xaxis.set_major_locator(monthloc)
                ax3.xaxis.set_major_formatter(monthfmt)
                # reduce fontsize of xtick labels 
                for label in ax3.xaxis.get_ticklabels():
                    label.set_fontsize(8)
                # reduce fontsize for y-axis labels
                for label in ax3.yaxis.get_ticklabels():
                    label.set_fontsize(8)
                for label in ax3d.yaxis.get_ticklabels():
                    label.set_fontsize(8)
                #put an expanded explanation in as a second subfigure
                #use siteinformation to add altitude, etc.
                ax2 = fig.add_subplot(212)  # make a handle (ax) to the axes object in the figure
                ax2.set_axis_off()
                ax2.text(0.5,0.70,'EXPLANATION',fontsize=10,ha='center',transform=ax2.transAxes)
                ax2.set_xlim([0.,1.])
                ax2.set_ylim([0.,1.])
                x1=0.2
                x2=0.25
                xscat=(x1+x2)/2.
                x3=0.30
                y1=0.60
                baseline=0.085
                fntsize=10
                ax2.plot([x1,x2],[y1,y1],'b-')
                ax2.text(x3,y1,"AVERAGE DAILY MAXIMUM DEPTH TO WATER",fontsize=fntsize,
                                 ha='left',va='center')
                y1=y1-baseline
                ax2.plot([x1,x2],[y1,y1],'r-')
                ax2.text(x3,y1,"AVERAGE DAILY MAXIMUM DEPTH TO WATER+STANDARD DEVIATION",fontsize=fntsize,
                                 ha='left',va='center')
                y1=y1-baseline
                ax2.plot([x1,x2],[y1,y1],'g-')
                ax2.text(x3,y1,"AVERAGE DAILY MAXIMUM DEPTH TO WATER-STANDARD DEVIATION",fontsize=fntsize,
                                 ha='left',va='center')
                y1=y1-baseline
                if currentyes==1:
                    ax2.plot([x1,x2],[y1,y1],'k-')
                    ax2.text(x3,y1,"DAILY MAXIMUM DEPTH TO WATER CURRENT WATER YEAR",fontsize=fntsize,
                                     ha='left',va='center')
                    y1=y1-baseline
                if tapeyes==1:
                    ax2.plot([xscat],[y1],'bx')
                    ax2.text(x3,y1,"TAPEDOWN DEPTH TO WATER CURRENT WATER YEAR",fontsize=fntsize,
                                     ha='left',va='center')
                    y1=y1-baseline
                
                y1=y1-baseline
                strg="USGS Site ID: {0:s}".format(cstation)
                ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                y1=y1-baseline
                strg="{0:s}".format(sitedat[cstation]['station_nm'])
                ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                y1=y1-baseline
                strg="Latitude, Longitude: {0:s}, {1:s}".format(sitedat[cstation]['dec_lat_va'],sitedat[cstation]['dec_long_va'])
                ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                y1=y1-baseline
                strg="Altitude of Land surface: {0:s} Feet {2:s}, Depth of Well: {1:s} Feet".format(sitedat[cstation]['alt_va'],
                                                                               sitedat[cstation]['well_depth_va'],
                                                                               sitedat[cstation]['alt_datum_cd'])
                ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                y1=y1-baseline
                #if local aquifer or national aquifer code is set, read in code dictionary
                #and use name instead of code in the explanation
                if codeflag==0:
                    if re.match('^\w',sitedat[cstation]['aqfr_cd'].strip()) or re.match('^\w',sitedat[cstation]['nat_aqfr_cd'].strip()):
                        code_file = os.path.join('..','NWIS_meta_data','codes.txt')
                        natl_code_file=os.path.join('..','NWIS_meta_data','national_codes.csv')
                        CODES=open(code_file,'r')
                        aqname=dict()
                        codeflag=1
                        for line in CODES:
                            vals=re.split('\t',line.strip())
                            aqname[vals[1]]=vals[2]
                        CODES.close()
                        NATLCODES=open(natl_code_file,'r')
                        NATLCODES.readline()
                        for line in NATLCODES:
                            vals=re.split(',',line.strip())
                            aqname[vals[1]]=vals[0]
                        NATLCODES.close()
                            
                if  sitedat[cstation]['aqfr_cd'] in aqname:
                    strg="Local Aquifer Code: {0:s}".format(aqname[sitedat[cstation]['aqfr_cd']].title())
                    ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                    y1=y1-baseline
                if sitedat[cstation]['nat_aqfr_cd'] in aqname:
                    strg="National Aquifer Code: {0:s}".format(aqname[sitedat[cstation]['nat_aqfr_cd']])
                    ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                    y1=y1-baseline
    ##            
                plt.savefig(os.path.join('figures',cstation+'_movingave' + plot_format),orientation='portrait')
                if disp_plot:
                    plt.show()

#--Function to plot hydrographs from an NWIS query using MATPLOTLIB
def NWIS_plotter(indat,gwdat,sitedat,plot_format,disp_plot=False):
    '''
    NWIS_plotter(indat,station_lookup,plot_format,disp_plot)
    
    A function to plot each hydrograph from an NWIS data query. The alternate name is used as the title
    Mike Fienen - 7/16/2012
    <mnfienen *at* usgs *dot* gov>

    Modified by HWR 8/16/2012 to plot both daily values and tapedown measurements, uses information
    from sitefile in generating plots
    
    INPUT:
    indat --> a dictionary of data returned from puller 'dv'
    gwdat --> dictionary of manual gw levels (tapedowns) from puller 'gwlevels'
    sitedat --> dictionary of sitefile information from puller 'site'
    plot_format --> format to save files of plots. can be '.png' or '.pdf'
    disp_plot --> a Boolean flag determining whether or not to display the plots onscreen
    
    OUTPUT:
    a set of files with plots of depth to water (reversing y axis) and named with station ID
    '''
    # check to see if there's already a figures subdirectory. If not, make one!
    if os.path.exists('figures')==False:
        os.mkdir('figures')

    codeflag=0
    for cstation in sitedat.keys():
        # tell the user what's going on
        print 'plotting Station ID: %s' %(cstation)
        #
        # first parse all the data for the current station
        #
        cname = sitedat[cstation]['station_nm'] # find the name from the stationID lookup dictionary
        dailyplot=0
        tapedownplot=0
        npcstation=np.array([cstation],dtype=str)
        if npcstation[0] in indat:
            DTW = indat[cstation]['DTW'] # get the depth to water for the current station
            prov_code = indat[cstation]['prov_code'] # get the provisional/approved codes
            datesdv = indat[cstation]['dates'] # pull the dates as datetime objects
            p_inds = np.nonzero(prov_code=='p')[0]
            a_inds = np.nonzero(prov_code=='a')[0]
            dailyplot=1
        if cstation in gwdat:    
            tapedwn = gwdat[cstation]['levels']
            datestp = gwdat[cstation]['dates']
            tapedownplot=1
        #
        # then set up the figure and plot it
        #
        if tapedownplot==0 and dailyplot==0:
            #no data for site in sitefile
            continue
        else:
            fig = plt.figure() # make a figure
            ax = fig.add_subplot(211)  # make a handle (ax) to the axes object in the figure
            #plt.hold = True  # hold the figure to allow multiple plotting events
            if dailyplot==1:
                if len(a_inds) > 0:
                    time_series_segments(datesdv,DTW,a_inds,'b-',10)
                if len(p_inds) > 0:
                    time_series_segments(datesdv,DTW,p_inds,'r-',10) 
                    
            if tapedownplot==1:
                plt.plot(datestp,tapedwn,'bx')

            #set the limits rounded to 5
            if dailyplot==1 and tapedownplot==1:
                ymin=min(np.min(DTW),np.min(tapedwn))
                ymax=max(np.max(DTW),np.max(tapedwn))
            elif dailyplot==1:
                ymin=np.min(DTW)
                ymax=np.max(DTW)
            elif tapedownplot==1:
                ymin=np.min(tapedwn)
                ymax=np.max(tapedwn)
            ymin=math.floor(ymin/5)*5.
            ymax=math.ceil(ymax/5)*5.
            ax.set_ylim([ymax,ymin])
            # set a second y-axis for altitude of water level
            # twinx() messed up x-tick labels, so need to do them
            # manually
            axd=ax.twinx()
            altmin=float(sitedat[cstation]['alt_va'])-ymax
            altmax=float(sitedat[cstation]['alt_va'])-ymin
            axd.set_ylim([altmin,altmax])
            # label the axes
            ax.set_ylabel('DEPTH TO WATER IN FEET',fontsize=10)
            ax.set_xlabel('DATE',fontsize=10)
            axd.set_ylabel('ALTITUDE OF WATER LEVEL IN FEET',fontsize=10)
            (xmin,xmax)=ax.get_xlim()
            delx=(xmax-xmin)/365
            yesmin=1
            if delx > 20:
                yearloc=mdates.YearLocator(2)
                minloc=mdates.YearLocator()
            elif delx > 30:
                yearloc=mdates.YearLocator(4)
                minloc=mdates.YearLocator()
            elif delx > 50:
                yearloc=mdates.YearLocator(5)
                minloc=mdates.YearLocator()
            elif delx > 70:
                yearloc=mdates.YearLocator(10)
                minloc=mdates.YearLocator()
            else:
                yearloc=mdates.YearLocator()
                yesmin=0
            yearfmt=mdates.DateFormatter('%Y')
            ax.xaxis.set_major_locator(yearloc)
            if yesmin==1:
                ax.xaxis.set_minor_locator(minloc)
            ax.xaxis.set_major_formatter(yearfmt)
            # rotate the xtick labels to avoid dates overlapping
            for label in ax.xaxis.get_ticklabels():
                label.set_rotation(30)
                label.set_fontsize(8)
            # reduce fontsize for y-axis labels
            for label in ax.yaxis.get_ticklabels():
                label.set_fontsize(8)
            for label in axd.yaxis.get_ticklabels():
                label.set_fontsize(8)
            #plt.gcf().subplots_adjust(bottom=0.1) # make a little room for the date labels

            #put an expanded explanation in as a second subfigure
            #use siteinformation to add altitude, etc.
            ax2 = fig.add_subplot(212)  # make a handle (ax) to the axes object in the figure
            ax2.set_axis_off()
            ax2.text(0.5,0.70,'EXPLANATION',fontsize=10,ha='center',transform=ax2.transAxes)
            ax2.set_xlim([0.,1.])
            ax2.set_ylim([0.,1.])
            x1=0.3
            x2=0.35
            xscat=(x1+x2)/2.
            x3=0.40
            y1=0.60
            baseline=0.085
            fntsize=10
            if dailyplot==1:
                if len(p_inds) > 0 and len(a_inds) > 0:
                    ax2.plot([x1,x2],[y1,y1],'b-')
                    ax2.text(x3,y1,"APPROVED DAILY DEPTH TO WATER",fontsize=fntsize,
                             ha='left',va='center')
                    y1=y1-baseline
                    ax2.plot([x1,x2],[y1,y1],'r-')
                    ax2.text(x3,y1,"PROVISIONAL DAILY DEPTH TO WATER",fontsize=fntsize,
                             ha='left',va='center')
                    y1=y1-baseline
                elif len(p_inds) > 0:
                    ax2.plot([x1,x2],[y1,y1],'r-')
                    ax2.text(x3,y1,"PROVISIONAL DAILY DEPTH TO WATER",fontsize=fntsize,
                             ha='left',va='center')
                    y1=y1-baseline
                elif len(a_inds) > 0:
                    ax2.plot([x1,x2],[y1,y1],'b-')
                    ax2.text(x3,y1,"APPROVED DAILY DEPTH TO WATER",fontsize=fntsize,
                             ha='left',va='center')
                    y1=y1-baseline
            if tapedownplot==1:
                ax2.plot([xscat],[y1],'bx')
                ax2.text(x3,y1,"TAPEDOWN DEPTH TO WATER",fontsize=fntsize,
                         ha='left',va='center')
                y1=y1-baseline

            y1=y1-baseline
            strg="USGS Site ID: {0:s}".format(cstation)
            ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
            y1=y1-baseline
            strg="{0:s}".format(sitedat[cstation]['station_nm'])
            ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
            y1=y1-baseline
            strg="Latitude, Longitude: {0:s}, {1:s}".format(sitedat[cstation]['dec_lat_va'],sitedat[cstation]['dec_long_va'])
            ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
            y1=y1-baseline
            strg="Altitude of Land surface: {0:s} Feet {2:s}, Depth of Well: {1:s} Feet".format(sitedat[cstation]['alt_va'],
                                                                               sitedat[cstation]['well_depth_va'],
                                                                               sitedat[cstation]['alt_datum_cd'])
            ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
            y1=y1-baseline
            #if local aquifer or national aquifer code is set, read in code dictionary
            #and use name instead of code in the explanation
            if codeflag==0:
                if re.match('^\w',sitedat[cstation]['aqfr_cd'].strip()) or re.match('^\w',sitedat[cstation]['nat_aqfr_cd'].strip()):
                    code_file = os.path.join('..','NWIS_meta_data','codes.txt')
                    natl_code_file=os.path.join('..','NWIS_meta_data','national_codes.csv')
                    CODES=open(code_file,'r')
                    aqname=dict()
                    codeflag=1
                    for line in CODES:
                        vals=re.split('\t',line.strip())
                        aqname[vals[1]]=vals[2]
                    CODES.close()
                    NATLCODES=open(natl_code_file,'r')
                    NATLCODES.readline()
                    for line in NATLCODES:
                        vals=re.split(',',line.strip())
                        aqname[vals[1]]=vals[0]
                    NATLCODES.close()
                        
            if  sitedat[cstation]['aqfr_cd'] in aqname:
                strg="Local Aquifer Code: {0:s}".format(aqname[sitedat[cstation]['aqfr_cd']].title())
                ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                y1=y1-baseline
            if sitedat[cstation]['nat_aqfr_cd'] in aqname:
                strg="National Aquifer Code: {0:s}".format(aqname[sitedat[cstation]['nat_aqfr_cd']])
                ax2.text(0.5,y1,strg,fontsize=fntsize,ha='center',va='center')
                y1=y1-baseline
##            
            plt.savefig(os.path.join('figures',cstation + plot_format),orientation='portrait')
            if disp_plot:
                plt.show()
        

if __name__ == '__main__':
    #set some values for testing and call the functions
    cStations=['454427084424002','434103083130301']
    puller_by_stations('site',cStations,'1950-10-01','2012-08-16','testsite.out')
#    (cState,cCounty)=['26','063']
#    puller_by_state_county('site',cState,cCounty,'1950-10-01','2010-09-30','testsiteSC.out','county')
    puller_by_stations('gwlevels',cStations,'1950-10-01','2012-08-16','tapedown.out')
#    puller_by_state_county('gwlevels',cState,cCounty,'1950-10-01','2010-09-30','tapedownSC.out','county')
    puller_by_stations('dv',cStations,'1950-10-01','2012-08-16','dv.out')
#    puller_by_state_county('dv',cState,cCounty,'1950-10-01','2010-09-30','dvSC.out','county')
# read files generated by puller
    infile='dv.out'
    (dvinfo,station_lookup)=NWIS_dv_reader(infile)
    infile='testsite.out'
    siteinfo=NWIS_sitefile_reader(infile)
    infile='tapedown.out'
    gwlevels=NWIS_gwlevels_reader(infile)
    NWIS_plotter(dvinfo,gwlevels,siteinfo,'.png',False)
    NWIS_plot_moving_average(dvinfo,gwlevels,siteinfo,'.png',False)
    
