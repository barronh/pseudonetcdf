#!/usr/bin/env python
import time
from datetime import datetime
import os
import argparse

import numpy as np
from netCDF4 import Dataset

def getdate(x):
    return datetime.strptime(x, '%Y-%m-%d').strftime('%Y%m%d')

parser = argparse.ArgumentParser("AQS Mart for pnceval")
parser.add_argument('-u', '--user', required = True, type = str, help = "User name for EPA's AQS Data Mart \"Query Air Data\" application at https://ofmext.epa.gov/AQDMRS/aqdmrs.html")
parser.add_argument('-p', '--password', required = True, type = str, help = "User name for EPA's AQS Data Mart \"Query Air Data\" application at https://ofmext.epa.gov/AQDMRS/aqdmrs.html")
parser.add_argument('-s', '--start-date', required = True, dest = 'bdate', type = getdate, help = 'Start date (inclusive) YYYYMMDD')
parser.add_argument('-e', '--end-date', required = True, dest = 'edate', type = getdate, help = 'Start date (inclusive) YYYYMMDD')
parser.add_argument('--param', type = str, default = '44201', nargs = '?', help = "Must exist as an AQS parameter")
parser.add_argument('--frmonly', type = str, choices = ['Y', 'N'], nargs = '?', default = 'Y', help = 'Federal register method only')
parser.add_argument('-d', '--download-output', dest = 'downout', type = str, nargs = '?', default = 'AQSREST_DATA.csv', help = 'Path for downloaded output.')
parser.add_argument('-o', '--output', type = str, dest = 'outpath', nargs = '?', default = 'AQSREST_DATA.nc', help = 'Path for output.')
parser.add_argument('-O', '--overwrite', dest = 'overwrite', default = False, action = 'store_true', help = 'Ovewrite if output already exists.')
parser.add_argument('GRIDCRO2D', help = 'CMAQ MCIP GRIDCRO2D file or any file that has LAT and LON variables')

args = parser.parse_args()


f = Dataset(args.GRIDCRO2D)
lon = f.variables['LON'][0, 0]
lat = f.variables['LAT'][0, 0]
args.minlon = llcrnrlon = lon[:, 0].max()
args.maxlon = urcrnrlon = lon[:, -1].min()
args.minlat = llcrnrlat = lat[0, :].max()
args.maxlat = urcrnrlat = lat[-1, :].min()
print args.minlon, args.maxlon
print args.minlat, args.maxlat
timezone = lon.mean() // 15.
if os.path.exists(args.outpath) and not args.overwrite:
    raise IOError('Path already exists: %s' % args.outpath)
    

def getrest():
    urlrequest = "https://ofmext.epa.gov/AQDMRS/ws/rawDataNotify?user=%(user)s&pw=%(password)s&format=DMCSV&param=%(param)s&bdate=%(bdate)s&edate=%(edate)s&minlat=%(minlat)s&maxlat=%(maxlat)s&minlon=%(minlon)s&maxlon=%(maxlon)s&dur=1&frmonly=%(frmonly)s" % dict(args._get_kwargs())


    import urllib2
    print urlrequest
    dataid = urllib2.urlopen(urlrequest).read()
    urlstatus = 'https://ofmext.epa.gov/AQDMRS/ws/status?id=' + dataid
    while True:
        status = urllib2.urlopen(urlstatus).read().strip()
        print status, dataid
        if status in ('Submitted', 'Processing'):
            time.sleep(3)
        elif status == 'Error':
            raise ValueError('AQS Mart did not correctly process your result.')
        elif status == 'Completed':
            break
        else:
            raise ValueError('AQS Mart returned unknown status (%s). Expected Submitted, Processing, Completed or Error.')

    urldata = 'https://ofmext.epa.gov/AQDMRS/ws/retrieve?id=' + dataid
    stmt = ('wget --no-check-certificate -O %s "%s"' % (args.downout, urldata))
    os.system(stmt)

#getrest()
#data = np.recfromtxt(args.downout, skip_footer = 2, delimiter = ",", names = True)
def hourly_parser(*args):
    import dateutil
    parse = np.vectorize(dateutil.parser.parse)
    out =  (args[0] + ' ' + args[1].astype('S2')+'Z').astype(np.datetime64)
    return out

import pandas
data = pandas.read_csv(args.downout, skip_footer = 2, parse_dates = [[12, 13]], date_parser = hourly_parser, converters = {'State Code': str, 'County Code': str, 'Site Num': str}, index_col = False, engine = 'python')

hourly = data.groupby(['Parameter Code', 'AQS Parameter Desc', 'Unit of Measure', 'Sample Duration', 'Date GMT_24 Hour GMT', 'State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean).sort(['Parameter Code', 'AQS Parameter Desc', 'Unit of Measure', 'Sample Duration', 'Date GMT_24 Hour GMT', 'State Code', 'County Code', 'Site Num'])
sites = data.groupby(['State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean)
nsites = len(sites)

outf = Dataset(args.outpath, 'w')
outf.createDimension('time', None)
outf.createDimension('LAY', 1)
sitelist = [row['State Code'] + row['County Code'] + row['Site Num'] for idx, row in sites.iterrows()]
outf.SITENAMES = ';'.join(sitelist)

outf.createDimension('points', len(sitelist))
outf.lonlatcoords = '/'.join(['%f,%f' % (row['Longitude'], row['Latitude']) for idx, row in sites.iterrows()])
last_time = start_time = hourly['Date GMT_24 Hour GMT'].values.min()
end_time = hourly['Date GMT_24 Hour GMT'].values.max()
ntimes = np.timedelta64(end_time - start_time, 'h').astype('i') + 1
alltimes = [np.timedelta64(i, 'h') + start_time for i in range(ntimes)]
temp = {}
lat = outf.createVariable('latitude', 'f', ('points',))
lat.units = 'degrees_north'
lat.standard_name = 'latitude'
lat[:] = sites['Latitude'].values
lon = outf.createVariable('longitude', 'f', ('points',))
lon.units = 'degrees_east'
lon.standard_name = 'longitude'
lon[:] = sites['Longitude'].values

for idx, row in hourly.iterrows():
    this_time = row['Date GMT_24 Hour GMT']
    val = row['Sample Measurement']
    unit = row['Unit of Measure'].strip()
    aqsid = row['State Code'] + row['County Code'] + row['Site Num']
    sidx = sitelist.index(aqsid)
    var_name = row['AQS Parameter Desc']
    if not var_name in outf.variables.keys():
        var = outf.createVariable(var_name, 'f', ('time', 'LAY', 'points'))
        var.units = unit
        var.standard_name = var_name
        temp[var_name] = np.zeros((ntimes, 1, nsites), dtype = 'f')
    tmpvar = temp[var_name]
    var = outf.variables[var_name]
    assert(var.units == unit)
    if last_time != this_time:
        last_time = this_time
        tidx = alltimes.index(np.datetime64(this_time))
    if tidx > ntimes:
        raise ValueError('Times (%d) exceed expected (%d)' % (tidx, ntimes))
    
    tmpvar[tidx, 0, sidx] = val

for varkey, tempvals in temp.iteritems():
    outf.variables[varkey][:] = tempvals
