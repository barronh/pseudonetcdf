from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from numpy import fromfile

import numpy as np

from datetime import datetime
import sys
import os

_geos_units = dict(ALBEDO = 'unitless',
CLDTOT = 'unitless',
EFLUX = 'W/m2',
EVAP = 'unitless',
FRLAKE = 'unitless',
FRLAND = 'unitless',
FRLANDIC = 'unitless',
FROCEAN = 'unitless',
GRN = 'unitless',
GWETROOT = 'unitless',
GWETTOP = 'unitless',
HFLUX = 'W/m2',
LAI = 'm2/m2',
LWGNET = 'W/m2',
LWTUP = 'W/m2',
PARDF = 'W/m2',
PARDR = 'W/m2',
PBLH = 'm',
PRECANV = 'kg/m2/s',
PRECCON = 'kg/m2/s',
PRECTOT = 'kg/m2/s',
PRECSNO = 'kg/m2/s',
QV2M = 'kg/kg',
SNOMAS = 'm',
SNODP = 'm',
SWGNET = 'W/m2',
T2M = 'K',
TROPP = 'hPa',
TSKIN = 'K',
U10M = 'm/s',
USTAR = 'm/s',
V10M = 'm/s',
Z0M = 'm',
DQRCON = 'kg/m2/s',
DQRLSC = 'kg/m2/s',
DTRAIN = 'kg/m2/s',
QL = 'unitless',
QI = 'unitless',
OPTDEPTH = 'unitless',
TAUCLI = 'unitless',
TAUCLW = 'unitless',
CLOUD = 'unitless',
PV = 'K*m2/kg/s',
OMEGA = 'Pa/s',
QV = 'kg/kg',
RH = 'unitless',
T = 'K',
U = 'm/s',
V = 'm/s',
DQIDTMST = 'kg/kg/s',
DQLDTMST = 'kg/kg/s',
DQVDTMST = 'kg/kg/s',
MOISTQ = 'kg/kg/s',
CMFMC = 'kg/m2/s', 	 	 	 	 
LWI = 'm/s',
PS = 'hPa',
SLP = 'hPa',
TO3 = 'DU',
TTO3 = 'DU')
class geos(PseudoNetCDFFile):
    def __init__(self, path, mode = 'r', reduced = True):
        from _vertcoord import geos_hyai, geos_hybi, geos_etam, geos_etai
        infile = file(path, 'r')
        infile.seek(0,0)
        fsize = os.path.getsize(path)
        lasttype = ''
        datatypes = [('name', '>i4, >S8, >i4')]
        names = []
        first = True
        while infile.tell() < fsize:
            ipad = fromfile(infile, dtype = '>i', count = 1)
            dtype = '>S8' if ipad == 8 else '>f'
            if dtype == '>f':
               elem = ipad / 4
               assert((ipad % 4) == 0)
            else:
               elem = 1
            if first:
                pass
            else:
                if elem == (nrow * ncol + 2):
                   shape = (nrow, ncol)
                elif elem == (nrow * ncol * nlay + 2):
                   shape = (nlay, nrow, ncol)
                elif elem == (nrow * ncol * nlay_stag + 2):
                   shape = (nlay_stag, nrow, ncol)
                
            if dtype == '>S8':
                data = fromfile(infile, dtype = dtype, count = elem)
            else:
                infile.seek(ipad, 1)
            #data = fromfile(infile, dtype = dtype, count = elem)
            #date, time = data[:2]
            #if dtype != '>c':
            #    data = data[2:].reshape(shape)
            epad = fromfile(infile, dtype = '>i', count = 1)
            if lasttype == '>S8' and name[:2] != 'G5':
                datatypes.append((name, '>i4, >S8, >i4, >i4, >i4, >i4, %s>f, >i4' % str(tuple(shape))))
            else:
                name = data[0].strip()
                if first:
                    first = False
                    if name[:2] == 'G5':
                        etak = 'GEOS-5-' + ('REDUCED' if reduced else 'NATIVE')
                        # Ap [hPa]
                        self.Ap = geos_hyai[etak]
                        # Bp [unitless]
                        self.Bp = geos_hybi[etak]
                        eta_m = geos_etam[etak]
                        eta_i = geos_etai[etak]
                        nlay = self.Ap.size - 1
                        nlay_stag = self.Ap.size
                    elif name[:2] == 'G4':
                        etak = 'GEOS-4-' + ('REDUCED' if reduced else 'NATIVE')
                        # Ap [hPa]
                        self.Ap = geos_hyai[etak]
                        # Bp [unitless]
                        self.Bp = geos_hybi[etak]
                        eta_m = geos_etam[etak]
                        eta_i = geos_etai[etak]
                        nlay = self.Ap.size - 1
                        nlay_stag = self.Ap.size
                    else:
                        raise ValueError('Unknown type %s' % name)

                    if name[3:5] == '22':
                        longitude_bounds = np.arange(-181.25, 180, 2.5).repeat(2, 0)[1:-1].reshape(-1, 2)
                        latitude_bounds = np.append(np.append(-90., np.arange(-89., 90., 2.)), 90.).repeat(2, 0)[1:-1].reshape(-1, 2)
                        nrow = 91
                        ncol = 144
                    elif name[3:5] == '45':
                        longitude_bounds = np.arange(-182.5, 180, 5).repeat(2, 0)[1:-1].reshape(-1, 2)
                        latitude_bounds = np.append(np.append(-90., np.arange(-88., 90., 4.)), 90.).repeat(2, 0)[1:-1].reshape(-1, 2)
                        nrow = 46
                        ncol = 72
                    else:
                        raise ValueError('Unknown type %s' % name)

                    longitude = longitude_bounds.mean(1)
                    latitude = latitude_bounds.mean(1)
                if name in names:
                    break
                names.append(name)
            lasttype = dtype
            assert(ipad == epad)
        if 'a3' in path:
            nsteps = 8
        elif 'a6' in path or 'i6' in path:
            nsteps = 4
        datatype = [datatypes[0], ('data', np.dtype(datatypes[1:]), nsteps)]
        data = np.memmap(path, dtype = np.dtype(datatype), mode = mode)
        d = self.createDimension('time', nsteps)
        d.setunlimited(True)
        self.createDimension('latitude', nrow)
        self.createDimension('longitude', ncol)
        self.createDimension('layer', nlay)
        self.createDimension('layer_stag', nlay_stag)
        self.createDimension('nv', 2)
        self.title = data['name']['f1'][0].strip()
        def getem(key):
            thisblock = data[0]['data'][key]
            thisdata = thisblock['f6']
            assert((thisblock['f3'] == thisblock['f7']).all())
            if len(thisdata.shape) == 3:
                dims = ('time', 'latitude', 'longitude')
            elif thisdata.shape[1] == nlay:
                dims = ('time', 'layer', 'latitude', 'longitude')
            else:
                dims = ('time', 'layer_stag', 'latitude', 'longitude')
            unit = _geos_units.get(key, '')    
            return PseudoNetCDFVariable(self, key, 'f', dims, values = thisdata, units = unit, long_name = key.ljust(16))
        self.variables = PseudoNetCDFVariables(keys = names[1:], func = getem)
        dates = data[0]['data'][names[1]]['f4']
        times = data[0]['data'][names[1]]['f5']
        years, months, days = dates // 10000, dates % 10000 // 100, dates % 100
        hours, minutes, seconds = times // 10000, times % 10000 // 100, times % 100
        self.createVariable('layer', 'f', ('layer',), values = eta_m[:], units = 'hybrid pressure-sigma', bounds = 'layer_bounds')
        self.createVariable('layer_bounds', 'f', ('layer','nv'), values = eta_i[:].repeat(2, 0)[1:-1].reshape(-1, 2), units = 'hybrid pressure-sigma')
                
        self.createVariable('time', 'f', ('time',), values = np.array([(datetime(y, m, d, H, M, S) - datetime(1985, 1, 1, 0, 0, 0)).total_seconds() / 3600. for y, m, d, H, M, S in zip(years, months, days, hours, minutes, seconds)]), units = 'hours since 1985-01-01 00:00:00 UTC', base_units = 'hours since 1985-01-01 00:00:00 UTC', long_name = 'time', standard_name = 'time')
        self.createVariable('longitude', 'f', ('longitude',), values = longitude, units = 'degrees_east', standard_name = 'longitude', bounds = 'longitude_bounds')
        self.createVariable('longitude_bounds', 'f', ('longitude','nv'), values = longitude_bounds, units = 'degrees_east', standard_name = 'longitude')
        self.createVariable('latitude', 'f', ('latitude',), values = latitude, units = 'degrees_north', standard_name = 'latitude', bounds = 'latitude_bounds')
        self.createVariable('latitude_bounds', 'f', ('latitude','nv'), values = latitude_bounds, units = 'degrees_north', standard_name = 'latitude')

import unittest
from warnings import warn
class TestMemmaps(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testGEOS(self):
        warn('Test not implemented')

if __name__ == '__main__':
    path = sys.argv[1]
    f = geos(path)
