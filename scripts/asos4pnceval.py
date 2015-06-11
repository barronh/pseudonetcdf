#!/usr/bin/env python
from scipy.constants import F2K
import os
import json
import datetime
import urllib2


# IEM quirk to have Iowa AWOS sites in its own labeled network
from shapely.geometry import Polygon, Point
import numpy as np
import argparse
from netCDF4 import Dataset
parser = argparse.ArgumentParser(prog = 'METAR', description = "Downloads met data from http://mesonet.agron.iastate.edu for a CMAQ domain")
parser.add_argument('--network', default = 'US__ASOS', help = 'Network follwing country_state_type conventions e.g., US_FL_ASOS')
parser.add_argument('GRIDCRO2D', help = 'CMAQ MCIP GRIDCRO2D file or any file that has LAT and LON variables')
parser.add_argument('--START_DATE', default = None, help = 'YYYY-MM-DD')
parser.add_argument('--END_DATE', default = None, help = 'YYYY-MM-DD')
args = parser.parse_args()

f = Dataset(args.GRIDCRO2D)
lon = f.variables['LON'][0, 0]
lat = f.variables['LAT'][0, 0]
llcrnrlon = lon[:, 0].max()
urcrnrlon = lon[:, -1].min()
llcrnrlat = lat[0, :].max()
urcrnrlat = lat[-1, :].min()

if args.START_DATE is None:
    startts = datetime.datetime.strptime(str(f.SDATE), '%Y%j')
else:
    startts = datetime.datetime.strptime(args.START_DATE, '%Y-%m-%d')
if args.START_DATE is None:
    endts = datetime.datetime.strptime(str(f.SDATE), '%Y%j')
else:
    endts = datetime.datetime.strptime(args.END_DATE, '%Y-%m-%d')

bounds = Polygon([(llcrnrlon, llcrnrlat), (llcrnrlon, urcrnrlat), (urcrnrlon, urcrnrlat), (urcrnrlon, llcrnrlat)])

# timestamps in UTC with latlon in tab delimited format to request data for

SERVICE = "http://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
SERVICE += "data=all&tz=Etc/UTC&format=tdf&latlon=yes&"

SERVICE += startts.strftime('year1=%Y&month1=%m&day1=%d&')
SERVICE += endts.strftime('year2=%Y&month2=%m&day2=%d&')

sites = []
sitelocs = {}
for network in args.network.split(','):
    # Get metadata
    uri = "http://mesonet.agron.iastate.edu/geojson/network.php?network=%s" % (
                                                                    network,)
    data = urllib2.urlopen(uri)
    jdict = json.load(data)

    nsites = [f['properties']['sid'] for f in jdict['features'] if bounds.intersects(Point(f['geometry']['coordinates']))]
    nlocs = dict([(f['properties']['sid'], f['geometry']['coordinates']) for f in jdict['features'] if bounds.intersects(Point(f['geometry']['coordinates']))])
    print nsites
    for site in nsites:
        sitelocs[site] = nlocs[site]
        uri = '%s&station=%s' % (SERVICE, site)
        print 'Network: %s Downloading: %s' % (network, site)
        outfn = '%s_%s_%s.txt' % (site, startts.strftime("%Y%m%d%H%M"),
                                  endts.strftime("%Y%m%d%H%M"))
        if os.path.exists(outfn):
            print 'Using cached'
            continue
        data = urllib2.urlopen(uri)
        out = open(outfn, 'w')
        out.write(data.read())
        out.close()
    sites += nsites
    
    
outf = Dataset('%s_obs_%s_%s.nc' % (os.path.basename(args.GRIDCRO2D).replace('.nc', ''), startts.strftime('%Y-%m-%d'), endts.strftime('%Y-%m-%d')), mode = 'w')
outf.createDimension('TSTEP', None)
outf.createDimension('VAR', 4)
outf.createDimension('DATE-TIME', 2)
outf.createDimension('LAY', 1)
outf.createDimension('points', len(sites))
outf.SDATE = int(startts.strftime('%Y%j'))
outf.STIME = 0
time_bounds = [(startts + datetime.timedelta(hours = i), startts + datetime.timedelta(hours = i + 1)) for i in range(int((endts - startts).total_seconds() // 3600) + 24)]

tflag = outf.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
nstep = len(time_bounds)
tflag[0:nstep, :, 0] = 0
tflag[0:nstep, :, 0] = np.array([int(t[0].strftime('%Y%j')) for t in time_bounds])[:, None].repeat(4, 1)
tflag[0:nstep, :, 1] = np.array([int(t[0].strftime('%H%M%S')) for t in time_bounds])[:, None].repeat(4, 1)

time = outf.createVariable('time', 'f', ('TSTEP',), fill_value = -999)
time[:] = tflag[:, 0, 0] + tflag[:, 0, 1] / 240000.

temp2 = outf.createVariable('TEMP2', 'f', ('TSTEP', 'LAY', 'points'), fill_value = -999)
temp2.units = 'K'.ljust(16)

wdir10 = outf.createVariable('WDIR10', 'f', ('TSTEP', 'LAY', 'points'), fill_value = -999)
wdir10.units = 'degrees from N'.ljust(16)

wspd10 = outf.createVariable('WSPD10', 'f', ('TSTEP', 'LAY', 'points'), fill_value = -999)
wspd10.units = 'M/S'.ljust(16)

rn = outf.createVariable('RN', 'f', ('TSTEP', 'LAY', 'points'), fill_value = -999)
rn.units = 'CM'.ljust(16)

lat = outf.createVariable('LAT', 'f', ('points'))
lat.units = 'degrees_north'.ljust(16)
lon = outf.createVariable('LON', 'f', ('points'))
lon.units = 'degrees_east'.ljust(16)

for si, site in enumerate(sites):
    outfn = '%s_%s_%s.txt' % (site, startts.strftime("%Y%m%d%H%M"),
                                  endts.strftime("%Y%m%d%H%M"))
    data = np.recfromtxt(outfn, missing_values = 'M', delimiter = '\t', names = True, comments = '#', skiprows = 5, usemask = True)
    lat[si] = sitelocs[site][1]
    lon[si] = sitelocs[site][0]
    if data.size == 0: continue
    times = np.array([datetime.datetime.strptime(d, '%Y-%m-%d %H:%M') for d in data['valid']])
    for ti, (sh, eh) in enumerate(time_bounds):
        tdata = data[(times > sh) & (times <=  eh)]
        if 'tmpf' in data.dtype.names:
            val = F2K(tdata['tmpf']).mean()
            if not np.ma.is_masked(val):
                temp2[ti, 0, si] = val
        if 'drct' in data.dtype.names:
            val = tdata['drct'].mean()
            if not np.ma.is_masked(val):
                wdir10[ti, 0, si] = val
        if 'sknt' in data.dtype.names:
            val = tdata['sknt'].mean() * 0.514444444 # convert knot to m/s
            if not np.ma.is_masked(val):
                wspd10[ti, 0, si] = val
        if 'p01i' in data.dtype.names:
            val = tdata['p01i'].mean() * 2.54 # convert inches to cm
            if not np.ma.is_masked(val):
                rn[ti, 0, si] = val
    outf.sync()
outf.close()
