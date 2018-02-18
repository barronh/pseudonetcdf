#!/usr/bin/env python
from __future__ import print_function

import os

from dateutil.parser import parse
from datetime import datetime, timedelta
import argparse

from netCDF4 import Dataset
import numpy as np
try:
    import pandas
except:
    raise ImportError('pncaqsraw4pnceval requires pandas; install pandas (e.g., pip install pandas)')
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
    """, formatter_class = argparse.RawDescriptionHelpFormatter)
parser.add_argument('-v', '--verbose', action = 'count', dest = 'verbose', default = 0)
parser.add_argument('--sampleval', default = None, help = 'Defaults to "Sample Measurement" for hourly and "Arithmetic Mean" for daily')
parser.add_argument('--timeresolution', choices = ['daily', 'hourly'], default = 'hourly', help = 'Defaults to hourly')
parser.add_argument('-s', '--start-date', required = True, dest = 'bdate', type = getbdate, help = 'Start date (inclusive) YYYY-MM-DD')
parser.add_argument('-e', '--end-date', required = True, dest = 'edate', type = getedate, help = 'End date (inclusive) YYYY-MM-DD')
parser.add_argument('-r', '--ref-date', default = datetime(1900, 1, 1), dest = 'rdate', type = getrdate, help = 'Reference date YYYYMMDD HH:MM:SS')
parser.add_argument('--param', type = str, default = '44201', nargs = '?', help = "Must exist as an AQS parameter")
spacegroup = parser.add_mutually_exclusive_group(required = False)
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

years = np.arange(args.bdate.year, args.edate.year + 1)
yearpaths = []
for year in years:
    yearpath = '%s_%s_%s.csv' % (args.timeresolution, args.param, year)
    filepath = '%s_%s_%s.zip' % (args.timeresolution, args.param, year)
    print('Downloading', yearpath)
    url = 'http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/%s' % (filepath,)
    if os.path.exists(yearpath):
        print('..Already have ' + yearpath)
    else:
        if os.path.exists(filepath):
            print('..Already have ' + filepath)
        else:
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
    yearpaths.append(yearpath)

from PseudoNetCDF.epafiles import aqsraw
from PseudoNetCDF.sci_var import stack_files
from PseudoNetCDF.pncgen import pncgen

infiles = []
for yearpath in yearpaths:
    infile = aqsraw(yearpath, timeresolution = args.timeresolution, bdate = args.bdate, edate = args.edate, rdate = args.rdate, wktpolygon = args.wktpolygon, sampleval = args.sampleval, verbose = args.verbose)
    infiles.append(infile)
    
outfile = stack_files(infiles, 'time')
persisted = pncgen(outfile, args.outpath)
persisted.sync()
print('Successful')
