#!/usr/bin/env python
from __future__ import print_function

import os

from dateutil.parser import parse
from datetime import datetime, timedelta
import argparse

from netCDF4 import Dataset
import numpy as np
import pandas
from shapely.wkt import loads
from shapely.geometry import Point, Polygon
from shapely.prepared import prep

def getbdate(x):
    return datetime.strptime(x, '%Y-%m-%d')

def getedate(x):
    return datetime.strptime(x + ' 23:59', '%Y-%m-%d %H:%M')

def getrdate(x):
    return datetime.strptime(x, '%Y-%m-%d %H:%M:%S')

parser = argparse.ArgumentParser(description = """Converts AQS Raw Hourly files for comparison with pncgen --extract files.

Example Workflow:
    $ %s --start-date=2006-08-01 --end-date=2012-08-01 --param=44201 GRIDCRO2D_Benchmark
    $ pncdump --header AQS_DATA_20060801-20060801.nc | grep lonlatcoords
    		:lonlatcoords = "-87.881412,30.498001/-85.802182,33.281261/..." ;
    $ pncgen -s LAY,0 --extract="-87.881412,30.498001/-85.802182,33.281261/..." CCTM_V5g_par_Linux2_x86_64gfort.ACONC.CMAQ-BENCHMARK_20060801 Benchmark_20060801-20060801.nc
    $ pnceval AQS_DATA_20060801-20060801.nc Benchmark_20060801-20060801.nc
    """)
parser.add_argument('--sampleval', default = None, help = 'Defaults to "Sample Measurement" for hourly and "Arithmetic Mean" for daily')
parser.add_argument('--timeresolution', choices = ['daily', 'hourly'], default = 'hourly', help = 'Defaults to hourly')
parser.add_argument('-s', '--start-date', required = True, dest = 'bdate', type = getbdate, help = 'Start date (inclusive) YYYY-MM-DD')
parser.add_argument('-e', '--end-date', required = True, dest = 'edate', type = getedate, help = 'Start date (inclusive) YYYY-MM-DD')
parser.add_argument('-r', '--ref-date', default = datetime(1900, 1, 1), dest = 'rdate', type = getrdate, help = 'Reference date YYYYMMDD HH:MM:SS')
parser.add_argument('--param', type = str, default = '44201', nargs = '?', help = "Must exist as an AQS parameter")
spacegroup = parser.add_mutually_exclusive_group(required = True)
spacegroup.add_argument('--gridcro2d', dest = 'GRIDCRO2D', help = 'CMAQ MCIP GRIDCRO2D file or any file that has LAT and LON variables')
spacegroup.add_argument('--wktpolygon')
parser.add_argument('-o', '--output', type = str, dest = 'outpath', nargs = '?', help = 'Path for output defaults to AQS_DATA_YYYYMMDD-YYYYMMDD.nc.')
parser.add_argument('-O', '--overwrite', dest = 'overwrite', default = False, action = 'store_true', help = 'Ovewrite if output already exists.')

args = parser.parse_args()
if args.sampleval is None:
    if args.timeresolution == 'hourly':
        args.sampleval = "Sample Measurement"
    elif args.timeresolution == 'daily':
        args.sampleval = "Arithmetic Mean"
    else:
        raise KeyError()

if args.outpath is None:
    args.outpath = 'AQS_DATA_%s-%s.nc' % (args.bdate.strftime('%Y%m%d'), args.edate.strftime('%Y%m%d'))

if os.path.exists(args.outpath) and not args.overwrite:
    raise IOError('Path already exists: %s' % args.outpath)

if not args.GRIDCRO2D is None:
    f = Dataset(args.GRIDCRO2D)
    lon = f.variables['LON'][0, 0]
    lat = f.variables['LAT'][0, 0]
    args.minlon = llcrnrlon = lon[:, 0].max()
    args.maxlon = urcrnrlon = lon[:, -1].min()
    args.minlat = llcrnrlat = lat[0, :].max()
    args.maxlat = urcrnrlat = lat[-1, :].min()
    ll = (args.minlon, args.minlat)
    lr = (args.maxlon, args.minlat)
    ur = (args.maxlon, args.maxlat)
    ul = (args.minlon, args.maxlat)
    args.wktpolygon = "POLYGON ((%s %s, %s %s, %s %s, %s %s, %s %s))" % (ll + lr + ur + ul + ll)

bounds = loads(args.wktpolygon)
pbounds = prep(bounds)
        

nseconds = {'hourly': 3600, 'daily': 3600*24}[args.timeresolution]
tunit = {'hourly': 'hours', 'daily': 'days'}[args.timeresolution]
ntimes = int((args.edate - args.bdate).total_seconds() / nseconds) + 1
alltimes = [timedelta(**{tunit: i}) + args.bdate for i in range(ntimes)]


# http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/download_files.html#Raw

def hourly_parser(*args):
    import dateutil
    parse = np.vectorize(dateutil.parser.parse)
    out =  (args[0] + ' ' + args[1].astype('S2')+'Z').astype(np.datetime64)
    return out

def hourly_parser(*args):
    import dateutil
    parse = np.vectorize(dateutil.parser.parse)
    out =  (args[0] + ' ' + args[1].astype('S2')+'Z').astype(np.datetime64)
    return out

years = np.arange(args.bdate.year, args.edate.year + 1)
hourlys = []
for year in years:
    yearpath = '%s_%s_%s.csv' % (args.timeresolution, args.param, year)
    filepath = '%s_%s_%s.zip' % (args.timeresolution, args.param, year)
    print('Downloading', yearpath)
    url = 'http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/%s' % (filepath,)
    if not os.path.exists(yearpath):
        if not os.path.exists(filepath):
            from urllib.request import urlopen, Request
            req = Request(url)
            ret = urlopen(req)
            data = ret.read()
            zipout = open(filepath, 'wb')
            zipout.write(data)
            zipout.flush()
            zipout.close()
        
        import zipfile
        zf = zipfile.ZipFile(filepath)
        zf.extract(yearpath)        

    print('Reading', yearpath)
    if args.timeresolution == 'hourly':
        parse_dates = [[11, 12]]
        date_key = 'Date GMT_Time GMT'
    else:
        parse_dates = [11]
        date_key = 'Date Local'
    data = pandas.read_csv(yearpath, index_col = False, converters = {'State Code': str, 'County Code': str, 'Site Num': str}, parse_dates = parse_dates)
    inspace = np.array([pbounds.contains(Point(plon, plat)) for plon, plat in zip(data['Longitude'].values, data['Latitude'].values)])
    intime = (data[date_key] >= args.bdate) & \
       (data[date_key] < args.edate)
    mask = inspace & intime
       

    hourly = data[mask].groupby(['Parameter Code', 'Parameter Name', 'Units of Measure', date_key, 'State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean).sort_values(by = ['Parameter Code', 'Parameter Name', 'Units of Measure', date_key, 'State Code', 'County Code', 'Site Num'])
    
    hourlys.append(hourly)

print('Concatenating files')
if len(hourlys) > 1:
    hourly = pandas.concat(hourlys)
else:
    hourly = hourlys[0]



sites = hourly.groupby(['State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean)
nsites = len(sites)

print('Creating output file')
outf = Dataset(args.outpath, 'w')
outf.createDimension('time', None)
outf.createDimension('LAY', 1)
sitelist = [row['State Code'] + row['County Code'] + row['Site Num'] for idx, row in sites.iterrows()]
outf.SITENAMES = ';'.join(sitelist)

outf.createDimension('points', len(sitelist))
outf.lonlatcoords = '/'.join(['%f,%f' % (row['Longitude'], row['Latitude']) for idx, row in sites.iterrows()])
last_time = start_time = hourly[date_key].values.min()
end_time = hourly[date_key].values.max()
temp = {}
lat = outf.createVariable('latitude', 'f', ('points',))
lat.units = 'degrees_north'
lat.standard_name = 'latitude'
lat[:] = sites['Latitude'].values
lon = outf.createVariable('longitude', 'f', ('points',))
lon.units = 'degrees_east'
lon.standard_name = 'longitude'
lon[:] = sites['Longitude'].values

print('Processing data rows')

for idx, row in hourly.iterrows():
    this_time = row[date_key]
    val = row[args.sampleval]
    unit = row['Units of Measure'].strip()
    aqsid = row['State Code'] + row['County Code'] + row['Site Num']
    sidx = sitelist.index(aqsid)
    var_name = row['Parameter Name']
    if not var_name in outf.variables.keys():
        var = outf.createVariable(var_name, 'f', ('time', 'LAY', 'points'), fill_value = -999)
        var.units = unit
        var.standard_name = var_name
        temp[var_name] = np.ma.masked_all((ntimes, 1, nsites), dtype = 'f')
        temp[var_name].set_fill_value(-999)
    tmpvar = temp[var_name]
    var = outf.variables[var_name]
    assert(var.units == unit)
    if last_time != this_time:
        last_time = this_time
        tidx = alltimes.index(this_time.to_datetime())
        print('\b\b\b\b\b\b %3d%%' % int(tidx / float(ntimes) * 100), end = '', flush = True)
    if tidx > ntimes:
        raise ValueError('Times (%d) exceed expected (%d)' % (tidx, ntimes))
    
    tmpvar[tidx, 0, sidx] = val
print('')
print('Writing to disk')
for varkey, tempvals in temp.items():
    outf.variables[varkey][:] = tempvals

time = outf.createVariable('time', 'f', ('time',))
time.units = args.rdate.strftime(tunit + ' since %F')
time.standard_name = 'time'
time[:] = [(t - args.rdate).total_seconds() / nseconds for t in alltimes]
outf.sync()
print('Successful')