#!/usr/bin/env python
import os
import sys
import numpy as np
import datetime
from argparse import ArgumentParser
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from io import BytesIO
import gzip
import tempfile
import getpass
from netCDF4 import MFDataset, Dataset
from shapely.wkt import loads
from shapely.geometry import Point, Polygon
from shapely.prepared import prep

from PseudoNetCDF.sci_var import slice_dim, getvarpnc, stack_files
from PseudoNetCDF.pncgen import pncgen

parser = ArgumentParser()
parser.add_argument('-v', '--verbose', default = 0, action = 'count', help = 'Show progress')
parser.add_argument('--checkflags', type = str, default = 'CSVKQkG', help = 'MADIS DD flags (default passes at least stage 1 or in Accept list)')
parser.add_argument('--dtmin', dest = 'dt_min', default = -1800., help = 'Seconds from nominal time negative')
parser.add_argument('--dtmax', dest = 'dt_max', default = 1800., help = 'Seconds from nominal time positive')
parser.add_argument('--no-humditity', dest = 'humidity', default = True, action = 'store_false', help = 'Disable calculation of specific humidity')
parser.add_argument('--username', default = None, help = 'MADIS username')
parser.add_argument('--password', default = None, help = 'MADIS password')
parser.add_argument('--cache', default = False, action = 'store_true', help = 'cache MADIS files')
parser.add_argument('--type', choices = ["sao"   , "metar"  , "maritime", "mesonet", "raob"  , "acarsProfiles"  , "profiler", "profiler", "hydro" ], help = 'MADIS data types')
spacegroup = parser.add_mutually_exclusive_group(required = False)
spacegroup.add_argument('--gridcro2d', dest = 'GRIDCRO2D', help = 'CMAQ MCIP GRIDCRO2D file or any file that has LAT and LON variables')
spacegroup.add_argument('--wktpolygon')
parser.add_argument('--START_DATE', type = lambda x: datetime.datetime.strptime(x, '%Y-%m-%d %H:%M'), help = 'YYYY-MM-DD HH:MM')
parser.add_argument('--END_DATE', type = lambda x: datetime.datetime.strptime(x, '%Y-%m-%d %H:%M'), help = 'YYYY-MM-DD HH:MM')
parser.add_argument('--TSTEP', default = '1 hours')
parser.add_argument('--hourlypath', help = 'Path for optional hourly output.')
parser.add_argument('rawpath', help = 'Outpath for timeObs output')

args = parser.parse_args()

if not args.GRIDCRO2D is None:
    gridcrof = Dataset(args.GRIDCRO2D)
    lon = gridcrof.variables['LON'][0, 0]
    lat = gridcrof.variables['LAT'][0, 0]
    llcrnrlon = lon[:, 0].max()
    urcrnrlon = lon[:, -1].min()
    llcrnrlat = lat[0, :].max()
    urcrnrlat = lat[-1, :].min()
    bounds = Polygon([(llcrnrlon, llcrnrlat), (llcrnrlon, urcrnrlat), (urcrnrlon, urcrnrlat), (urcrnrlon, llcrnrlat)])
    prep_bounds = prep(bounds)
elif not args.wktpolygon is None:
    bounds = loads(args.wktpolygon)
    prep_bounds = prep(bounds)
else:
    bounds = None
    prep_bounds = None

level1 = dict(zip(("sao"   , "metar"  , "maritime", "mesonet", "raob"  , "acarsProfiles"  , "profiler", "profiler", "hydro"), ("point" , "point"  , "point"   , "LDAD"   , "point" , "point"          ,  "point"  , "LDAD"    , "LDAD")))[args.type]
level3 = dict(zip(("sao"   , "metar"  , "maritime", "mesonet", "raob"  , "acarsProfiles"  , "profiler", "profiler", "hydro"), ("netcdf", "netcdf" , "netcdf"  , "netCDF" , "netcdf", "netcdf"         ,  "netcdf" , "netCDF"  , "netCDF")))[args.type]

authvalues = dict(username = args.username, password = args.password)
auth = urlencode(authvalues)

template = 'ftp://pftp.madis-data.noaa.gov/archive/%%Y/%%m/%%d/%(level1)s/%(level2)s/%(level3)s/%%Y%%m%%d_%%H%%M.gz' % dict(level2 = args.type, level3 = level3, level1 = level1)

times = [args.START_DATE]
now = times[0]
inc = datetime.timedelta(hours = 1)
while now < args.END_DATE:
    now = now + inc
    times.append(now)



timefiles = {}
infiles = []
for time in times:
    if args.verbose > 0:
        print(time, 'start')
    url = time.strftime(template)
    cachedpath = time.strftime(os.path.join('madis', 'archive', '%Y', '%m', '%d', level1, args.type, level3, os.path.basename(url)))
    if not args.cache or not os.path.exists(cachedpath):
        if args.username is None:
            args.username = input('MADIS username:')

        if args.password is None:
            args.password = getpass.getpass(prompt='MADIS password: ', stream=None)
        if args.verbose > 0:
            print('downloading: ' + url)
        req = Request(url, data = auth)
        try:
            ftpf = urlopen(req)
        except:
            print('Failed to retrieve ' + url + '; continuing', file = sys.stderr)
            continue
        # create BytesIO temporary file
        if args.cache:
            dirpath = os.path.dirname(cachedpath)
            if not os.path.exists(dirpath):
                os.makedirs(dirpath)
            compressedFile = open(cachedpath, 'w+b')
        else:
            compressedFile = BytesIO()
        compressedFile.write(ftpf.read())
        compressedFile.flush()
    else:
        if args.verbose > 0: print('using cache ' + cachedpath)
        compressedFile = open(cachedpath, 'r+b')
    compressedFile.seek(0)
    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')
    diskf = tempfile.NamedTemporaryFile('w+b')  
    infiles.append(diskf)  
    diskf.write(decompressedFile.read())
    diskf.flush()
    if args.verbose > 0:
        print(time, 'end')


ncff = MFDataset([inf.name for inf in infiles], 'r')
lat = ncff.variables['latitude'][:]
lon = ncff.variables['longitude'][:]
points = zip(lon, lat)
if prep_bounds is None:
    found_point_ids = range(len(points))
else:
    found_point_ids = []
    for pi, point in enumerate(points):
        isin = prep_bounds.contains(Point(*point))
        if isin:
            found_point_ids.append(pi)
            if args.verbose > 1:
                print(point, isin)
        elif args.verbose > 2:
                print(point, isin)

varkeys = ['temperature', 'windDir', 'windSpeed', 'dewpoint', 'altimeter']
vardds = [k + 'DD' for k in varkeys]

if args.verbose > 1:
    print('Subset variables')

getvarkeys = varkeys + vardds + ['stationName', 'timeObs', 'timeNominal', 'elevation', 'latitude', 'longitude']

if args.verbose > 1:
    print('Slicing files')

from PseudoNetCDF.sci_var import Pseudo2NetCDF, PseudoNetCDFFile
p2p = Pseudo2NetCDF(verbose = 0)
outfile = PseudoNetCDFFile()
p2p.addDimensions(ncff, outfile)
outfile.createDimension('recNum', len(found_point_ids))
p2p.addGlobalProperties(ncff, outfile)

for vark in getvarkeys:
    p2p.addVariable(ncff, outfile, vark, data = False)

for vark in getvarkeys:
    invar = ncff.variables[vark]
    outvar = outfile.variables[vark]
    recid = list(invar.dimensions).index('recNum')
    outvar[:] = invar[:].take(found_point_ids, recid)

if args.humidity:
    varkeys.append('specificHumidity')

if args.verbose > 1:
    print('Making output files')

maxlen = len(outfile.dimensions['maxStaNamLen'])
stationNames = outfile.variables['stationName'][:].view('S%d' % maxlen).squeeze()
timeObs = outfile.variables['timeObs']
timeNominal = outfile.variables['timeNominal']

ustations = np.unique(stationNames)
utimes = np.sort(np.unique(timeObs[:]))
unomtimes = np.sort(np.unique(timeNominal[:]))

obstimefile = Dataset(args.rawpath, 'w', format = 'NETCDF4_CLASSIC')
obstimefile.createDimension('time', utimes.size)
obstimefile.createDimension('timeNominal', unomtimes.size)
obstimefile.createDimension('stations', ustations.size)
obstimefile.createDimension('maxStaNamLen', maxlen)

inelev = outfile.variables['elevation']
inlat = outfile.variables['latitude']
inlon = outfile.variables['longitude']

outstat = obstimefile.createVariable('stationNames', 'S1', ('stations', 'maxStaNamLen'))
outstat[:] = ustations.view('S1')

outlat = obstimefile.createVariable('latitude', 'f', ('stations',))
for pk in inlat.ncattrs():
    setattr(outlat, pk, getattr(inlat, pk))

outlon = obstimefile.createVariable('longitude', 'f', ('stations',))
for pk in inlon.ncattrs():
    setattr(outlon, pk, getattr(inlon, pk))

outelev = obstimefile.createVariable('elevation', 'f', ('stations',))
for pk in inelev.ncattrs():
    setattr(outelev, pk, getattr(inelev, pk))

if args.humidity:
    from PseudoNetCDF.derived.met import wmr_ptd
    outfile.variables['specificHumidity'] = wmr_ptd(outfile.variables['altimeter']/100., outfile.variables['dewpoint']-273.15, gkg = True)
    outfile.variables['specificHumidityDD'] = outfile.variables['dewpointDD']

station_idx = (stationNames[:, None] == ustations[None, :]).argmax(1)
time_idx = (timeObs[:][:, None] == utimes[None, :]).argmax(1)

time = obstimefile.createVariable('time', 'd', ('time',))
for pk in timeObs.ncattrs():
    setattr(time, pk, getattr(timeObs, pk))
time.units = time.units.replace('.0', '')

nomtime = obstimefile.createVariable('timeNominal', 'd', ('timeNominal',))
for pk in timeNominal.ncattrs():
    setattr(nomtime, pk, getattr(timeNominal, pk))
nomtime.units = nomtime.units.replace('.0', '')

for vark in varkeys:
    invar = outfile.variables[vark]
    outvar = obstimefile.createVariable(vark, invar.dtype, ('time', 'stations'), fill_value = -999)
    for pk in invar.ncattrs():
        setattr(outvar, pk, getattr(invar, pk))

tmpvals = np.zeros_like(outlat)
tmpvals[station_idx] = inlat[:]
outlat[:] = tmpvals[:]
tmpvals = np.zeros_like(outlon)
tmpvals[station_idx] = inlon[:]
outlon[:] = tmpvals[:]
tmpvals = np.zeros_like(outelev)
tmpvals[station_idx] = inelev[:]
outelev[:] = tmpvals[:]
time[:] = utimes[:]
nomtime[:] = unomtimes

for vark in varkeys:
    invar = outfile.variables[vark]
    outvar = obstimefile.variables[vark]
    indd = np.char.decode(outfile.variables[vark + 'DD'][:], 'ascii')
    inval = invar[:]
    tmpvals = np.zeros_like(outvar[:])
    passes = np.sum([indd[:] == flag for flag in args.checkflags], axis = 0) > 0
    fails = ~passes
    tmpvals[time_idx, station_idx] = np.ma.masked_where(fails, invar[:])
    outvar[:] = tmpvals

if args.verbose > 1:
    print('Saving hourly data')


if not args.hourlypath is None:
    dt_min = args.dt_min
    dt_max = args.dt_max

    from PseudoNetCDF.sci_var import Pseudo2NetCDF, PseudoNetCDFFile
    p2p = Pseudo2NetCDF(verbose = False)
    nomtimefile = PseudoNetCDFFile()
    p2p.addDimensions(obstimefile, nomtimefile)
    nomtimefile.createDimension('time', unomtimes.size)
    p2p.addGlobalProperties(obstimefile, nomtimefile)

    for vark in varkeys:
        p2p.addVariable(obstimefile, nomtimefile, vark, data = False)
        outvar = nomtimefile.variables[vark]
        var = obstimefile.variables[vark]
        timem = (var[:]/var[:])*time[:][:, None]
        dt = np.ma.MaskedArray(unomtimes[:, None, None]) - timem[None,:, :]
        nomtime_idx, points_idx = np.indices(dt[:, 0].shape)
        time_idx = np.abs(dt).argmin(1)
        outvals = var[:][time_idx, points_idx]
        min_dt = dt[nomtime_idx, time_idx, points_idx]
        dt_mask = (min_dt < dt_min) | (min_dt > dt_max)
        outvals = np.ma.masked_where(dt_mask, outvals)
        outvar[:] = outvals[:]

    pncgen(nomtimefile, args.hourlypath, verbose = 0)
# :DD_long_name = "QC data descriptor model:  QC summary values" ;
# :DD_reference = "AWIPS Technique Specification Package (TSP) 88-21-R2" ;
# :DD_values = "Z,C,S,V,X,Q,K,k,G, or B" ;
# :DD_value_Z = "No QC applied" ;
# :DD_value_C = "Passed QC stage 1" ;
# :DD_value_S = "Passed QC stages 1 and 2" ;
# :DD_value_V = "Passed QC stages 1, 2 and 3" ;
# :DD_value_X = "Failed QC stage 1" ;
# :DD_value_Q = "Passed QC stage 1, but failed stages 2 or 3 " ;
# :DD_value_K = "Passed QC stages 1, 2, 3, and 4" ;
# :DD_value_k = "Passed QC stage 1,2, and 3, failed stage 4 " ;
# :DD_value_G = "Included in accept list" ;
# :DD_value_B = "Included in reject list" ;

