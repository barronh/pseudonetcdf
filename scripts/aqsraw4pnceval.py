#!/usr/bin/env python
import os

from dateutil.parser import parse
from datetime import datetime, timedelta
import argparse

from netCDF4 import Dataset
import numpy as np
import pandas

def getbdate(x):
    return datetime.strptime(x, '%Y-%m-%d')

def getedate(x):
    return datetime.strptime(x + ' 23:59', '%Y-%m-%d %H:%M')

def getrdate(x):
    return datetime.strptime(x, '%Y-%m-%d %H:%M')

parser = argparse.ArgumentParser(description = """Converts AQS Raw Hourly files for comparison with pncgen --extract files.

Example Workflow:
    $ %s --start-date=2006-08-01 --end-date=2012-08-01 --param=44201 GRIDCRO2D_Benchmark
    $ pncdump --header AQS_DATA_20060801-20060801.nc | grep lonlatcoords
    		:lonlatcoords = "-87.881412,30.498001/-85.802182,33.281261/..." ;
    $ pncgen -s LAY,0 --extract="-87.881412,30.498001/-85.802182,33.281261/..." CCTM_V5g_par_Linux2_x86_64gfort.ACONC.CMAQ-BENCHMARK_20060801 Benchmark_20060801-20060801.nc
    $ pnceval AQS_DATA_20060801-20060801.nc Benchmark_20060801-20060801.nc
    """)
parser.add_argument('-s', '--start-date', required = True, dest = 'bdate', type = getbdate, help = 'Start date (inclusive) YYYYMMDD')
parser.add_argument('-e', '--end-date', required = True, dest = 'edate', type = getedate, help = 'Start date (inclusive) YYYYMMDD')
parser.add_argument('-r', '--ref-date', default = datetime(1900, 1, 1), dest = 'rdate', type = getrdate, help = 'Reference date YYYYMMDD HH:MM:DD')
parser.add_argument('--param', type = str, default = '44201', nargs = '?', help = "Must exist as an AQS parameter")
parser.add_argument('GRIDCRO2D', help = 'CMAQ MCIP GRIDCRO2D file or any file that has LAT and LON variables; used to identify domain of interest.')
parser.add_argument('-o', '--output', type = str, dest = 'outpath', nargs = '?', help = 'Path for output defaults to AQS_DATA_YYYYMMDD-YYYYMMDD.nc.')
parser.add_argument('-O', '--overwrite', dest = 'overwrite', default = False, action = 'store_true', help = 'Ovewrite if output already exists.')

args = parser.parse_args()

if args.outpath is None:
    args.outpath = 'AQS_DATA_%s-%s.nc' % (args.bdate.strftime('%Y%m%d'), args.edate.strftime('%Y%m%d'))

if os.path.exists(args.outpath) and not args.overwrite:
    raise IOError('Path already exists: %s' % args.outpath)

f = Dataset(args.GRIDCRO2D)
lon = f.variables['LON'][0, 0]
lat = f.variables['LAT'][0, 0]
args.minlon = llcrnrlon = lon[:, 0].max()
args.maxlon = urcrnrlon = lon[:, -1].min()
args.minlat = llcrnrlat = lat[0, :].max()
args.maxlat = urcrnrlat = lat[-1, :].min()
ntimes = int((args.edate - args.bdate).total_seconds() / 3600) + 1
alltimes = [timedelta(hours = i) + args.bdate for i in range(ntimes)]

print args.minlon, args.maxlon
print args.minlat, args.maxlat

# http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/download_files.html#Raw

def hourly_parser(*args):
    import dateutil
    parse = np.vectorize(dateutil.parser.parse)
    out =  (args[0] + ' ' + args[1].astype('S2')+'Z').astype(np.datetime64)
    return out

years = np.arange(args.bdate.year, args.edate.year + 1)
hourlys = []
for year in years:
    yearpath = 'hourly_%s_%s.csv' % (args.param, year)
    print 'Downloading', yearpath
    if not os.path.exists(yearpath):
        getcmd = 'wget --continue http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/hourly_%s_%s.zip' % (args.param, year)
        os.system(getcmd)
        os.system('unzip hourly_%s_%s.zip' % (args.param, year))
    print 'Reading', yearpath
    data = pandas.read_csv(yearpath, index_col = False, converters = {'State Code': str, 'County Code': str, 'Site Num': str}, parse_dates = [[11, 12]], date_parser = hourly_parser)

    mask = (data['Latitude'].values >= args.minlat) & (data['Latitude'].values <= args.maxlat) & \
       (data['Longitude'].values >= args.minlon) & (data['Longitude'].values <= args.maxlon) & \
       (data['Date GMT_Time GMT'] >= args.bdate) & \
       (data['Date GMT_Time GMT'] < args.edate)
       

    hourly = data[mask].groupby(['Parameter Code', 'Parameter Name', 'Units of Measure', 'Date GMT_Time GMT', 'State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean).sort(['Parameter Code', 'Parameter Name', 'Units of Measure', 'Date GMT_Time GMT', 'State Code', 'County Code', 'Site Num'])
    
    hourlys.append(hourly)

print 'Concatenating files'
if len(hourlys) > 1:
    hourly = pandas.concat(hourlys)
else:
    hourly = hourlys[0]



sites = hourly.groupby(['State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean)
nsites = len(sites)

print 'Creating output file'
outf = Dataset(args.outpath, 'w')
outf.createDimension('time', None)
outf.createDimension('LAY', 1)
sitelist = [row['State Code'] + row['County Code'] + row['Site Num'] for idx, row in sites.iterrows()]
outf.SITENAMES = ';'.join(sitelist)

outf.createDimension('points', len(sitelist))
outf.lonlatcoords = '/'.join(['%f,%f' % (row['Longitude'], row['Latitude']) for idx, row in sites.iterrows()])
last_time = start_time = hourly['Date GMT_Time GMT'].values.min()
end_time = hourly['Date GMT_Time GMT'].values.max()
temp = {}
lat = outf.createVariable('latitude', 'f', ('points',))
lat.units = 'degrees_north'
lat.standard_name = 'latitude'
lat[:] = sites['Latitude'].values
lon = outf.createVariable('longitude', 'f', ('points',))
lon.units = 'degrees_east'
lon.standard_name = 'longitude'
lon[:] = sites['Longitude'].values

print 'Processing data rows'

for idx, row in hourly.iterrows():
    this_time = row['Date GMT_Time GMT']
    val = row['Sample Measurement']
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
        print '\b\b\b\b\b\b %3d%%' % int(tidx / float(ntimes) * 100),
    if tidx > ntimes:
        raise ValueError('Times (%d) exceed expected (%d)' % (tidx, ntimes))
    
    tmpvar[tidx, 0, sidx] = val
print
print 'Writing to disk'
for varkey, tempvals in temp.iteritems():
    outf.variables[varkey][:] = tempvals

time = outf.createVariable('time', 'f', ('time',))
time.units = args.rdate.strftime('hours since %F UTC')
time.standard_name = 'time'
time[:] = [(t - args.rdate).total_seconds() / 3600 for t in alltimes]
outf.sync()
print 'Successful'