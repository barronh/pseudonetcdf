from __future__ import print_function
__all__ = ['cdtoms', 'tomsl3']


from PseudoNetCDF import PseudoNetCDFFile
from PseudoNetCDF.pncwarn import warn
from re import compile
import numpy as np
from datetime import datetime


dayre = compile(' Day:\s+(?P<jday>\d+) (?P<daystring>.{12})\s+EP/' +
                'TOMS CORRECTED OZONE GEN:\d+\.\d+\sV\d ALECT:\s+\d+' +
                ':\d+ [AP]M ')
edgere = compile(
    '.+(?P<bins>\d+)\s+bins.*' +
    '\s+[-+]?(?P<start>\d+\.\d+)\s*(?P<startdir>[WENS]).*' +
    '\s+[-+]?(?P<end>\d+\.\d+)\s*(?P<enddir>[WENS])\s+.*' +
    '\((?P<step>\d+.\d+)\s*degree steps.+')


def _groupdict(reo, line):
    result = reo.match(line)
    if result is None:
        return {}
    else:
        return result.groupdict()


def cdtoms(path, outfile=None):
    if outfile is None:
        outfile = PseudoNetCDFFile()
    inlines = open(path, 'r').readlines()
    dayline = inlines[0]
    daygrp = _groupdict(dayre, dayline)

    sdate = 0
    if 'daystring' in daygrp:
        date = datetime.strptime(daygrp['daystring'], '%b %d, %Y')
        rdate = datetime(1970, 1, 1)
        sdate = (date - rdate).total_seconds()
    else:
        import pandas as pd
        dayparts = dayline.split(' ')
        for i in [3, 2, 1]:
            daystr = ' '.join(dayparts[:i])
            try:
                date = pd.to_datetime(daystr, box=False)
                break
            except Exception as e:
                print(e)
        else:
            date = np.datetime64('1970-01-01')
        rdate = np.datetime64('1970-01-01')
        sdate = (date - rdate).astype('d') / 1e9

    longrp = _groupdict(edgere, inlines[1])
    latgrp = _groupdict(edgere, inlines[2])

    for propdict in [daygrp, longrp, latgrp]:
        for k, v in propdict.items():
            try:
                v = eval(v)
            except Exception:
                pass
            setattr(outfile, k, v)
    outfile.HISTORY = ''.join(inlines[:3])

    blat = latgrp.get('start', '59.5')
    bsn = latgrp.get('startdir', 'S')
    elat = latgrp.get('end', '59.5')
    esn = latgrp.get('enddir', 'N')
    latstep = float(latgrp.get('step', '1'))
    blat = {'N': 1, 'S': -1}[bsn] * float(blat)
    elat = {'N': 1, 'S': -1}[esn] * float(elat)

    blon = longrp.get('start', '179.375')
    bwe = longrp.get('startdir', 'W')
    elon = longrp.get('end', '179.375')
    ewe = longrp.get('enddir', 'E')
    lonstep = float(longrp.get('step', '1.25'))
    blon = {'E': 1, 'W': -1}[bwe] * float(blon)
    elon = {'E': 1, 'W': -1}[ewe] * float(elon)

    datalines = inlines[3:]
    lats = []
    for i, line in enumerate(datalines):
        if 'lat' not in line:
            datalines[i] = line[1:-1].rstrip()
        else:
            data, lat = line.split('lat =')
            datalines[i] = data[1:-1].rstrip()
            lats.append(lat.strip())

    nlats = len(lats)
    datablock = ''.join(datalines).replace(' ', '0')
    nlons = len(datablock) // 3 // nlats
    outfile.createDimension('time', 1)
    outfile.createDimension('latitude', nlats)
    outfile.createDimension('longitude', nlons)
    outfile.createDimension('nv', 2)

    var = outfile.createVariable('time', 'f', ('time',))
    var.units = 'seconds since 1970-01-01 00:00:00+0000'
    var[:] = sdate

    var = outfile.createVariable('latitude', 'f', ('latitude',))
    var.units = 'degrees N'
    var[:] = np.arange(blat, elat + latstep, latstep)
    lat = var
    linelat = np.array(lats, dtype='f')
    if not (lat[:] == linelat).all():
        warn('Header metadata does not match lats')
        lat[:] = linelat

    var = outfile.createVariable('latitude_bounds', 'f', ('latitude', 'nv'))
    var.units = 'degrees N'
    var[:, 0] = lat - latstep / 2.
    var[:, 1] = lat + latstep / 2.

    var = outfile.createVariable('longitude', 'f', ('longitude',))
    var.units = 'degrees E'
    lon = var[:] = np.arange(blon, elon + lonstep, lonstep)

    var = outfile.createVariable('longitude_bounds', 'f', ('longitude', 'nv'))
    var.units = 'degrees E'
    var[:, 0] = lon - lonstep / 2.
    var[:, 1] = lon + lonstep / 2.

    var = outfile.createVariable(
        'ozone', 'f', ('latitude', 'longitude'),
        missing_value=999
    )
    var.units = 'matm-cm'
    var.long_name = var.var_desc = 'ozone'.ljust(16)
    var[:] = np.ma.masked_values(
        np.array(
            [i for i in datablock],
            dtype='S1'
        ).view('S3').astype('i'),
        var.missing_value
    ).reshape(nlats, nlons)

    return outfile


class tomsl3(PseudoNetCDFFile):
    def __init__(self, path):
        cdtoms(path, self)


if __name__ == '__main__':
    from PseudoNetCDF.pncdump import pncdump
    pfile = cdtoms('test.txt')
    pncdump(pfile)
    print(pfile.variables['ozone'][2, 3])
