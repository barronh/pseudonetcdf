from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from numpy import fromfile

import numpy as np

from datetime import datetime
import sys
import os
g4_eta_z = np.array([[  1.00000000e+00,  -6.00000000e-03],
       [  9.92555000e-01,   5.80000000e-02],
       [  9.85110000e-01,   1.22000000e-01],
       [  9.70469000e-01,   2.49000000e-01],
       [  9.55828000e-01,   3.78000000e-01],
       [  9.29330000e-01,   6.15000000e-01],
       [  9.02831000e-01,   8.57000000e-01],
       [  8.66492000e-01,   1.19700000e+00],
       [  8.30152000e-01,   1.54900000e+00],
       [  7.86601000e-01,   1.98600000e+00],
       [  7.43050000e-01,   2.44300000e+00],
       [  6.95218000e-01,   2.96900000e+00],
       [  6.47387000e-01,   3.52300000e+00],
       [  5.98479000e-01,   4.12500000e+00],
       [  5.49571000e-01,   4.76800000e+00],
       [  5.08017000e-01,   5.35100000e+00],
       [  4.66464000e-01,   5.97500000e+00],
       [  4.31179000e-01,   6.54100000e+00],
       [  3.95895000e-01,   7.14800000e+00],
       [  3.65938000e-01,   7.70000000e+00],
       [  3.35980000e-01,   8.29100000e+00],
       [  3.10561000e-01,   8.83000000e+00],
       [  2.85142000e-01,   9.40900000e+00],
       [  2.63587000e-01,   9.93600000e+00],
       [  2.42032000e-01,   1.05040000e+01],
       [  2.23772000e-01,   1.10210000e+01],
       [  2.05513000e-01,   1.15780000e+01],
       [  1.90061000e-01,   1.20860000e+01],
       [  1.74608000e-01,   1.26330000e+01],
       [  1.61513000e-01,   1.31340000e+01],
       [  1.48418000e-01,   1.36740000e+01],
       [  1.37287000e-01,   1.41700000e+01],
       [  1.26157000e-01,   1.47060000e+01],
       [  1.16695000e-01,   1.51980000e+01],
       [  1.07234000e-01,   1.57310000e+01],
       [  9.91910000e-02,   1.62220000e+01],
       [  9.11490000e-02,   1.67530000e+01],
       [  8.43130000e-02,   1.72430000e+01],
       [  7.74770000e-02,   1.77730000e+01],
       [  7.16000000e-02,   1.82690000e+01],
       [  6.57230000e-02,   1.88070000e+01],
       [  6.06820000e-02,   1.93090000e+01],
       [  5.56410000e-02,   1.98550000e+01],
       [  5.13260000e-02,   2.03640000e+01],
       [  4.70120000e-02,   2.09200000e+01],
       [  4.33260000e-02,   2.14380000e+01],
       [  3.96410000e-02,   2.20040000e+01],
       [  3.64990000e-02,   2.25310000e+01],
       [  3.33580000e-02,   2.31080000e+01],
       [  3.06730000e-02,   2.36480000e+01],
       [  2.79870000e-02,   2.42400000e+01],
       [  2.56990000e-02,   2.47940000e+01],
       [  2.34100000e-02,   2.54020000e+01],
       [  2.14670000e-02,   2.59710000e+01],
       [  1.95230000e-02,   2.65960000e+01],
       [  1.78780000e-02,   2.71800000e+01],
       [  1.62320000e-02,   2.78240000e+01],
       [  1.48440000e-02,   2.84230000e+01],
       [  1.34550000e-02,   2.90850000e+01],
       [  1.22870000e-02,   2.97010000e+01],
       [  1.11200000e-02,   3.03830000e+01],
       [  1.01410000e-02,   3.10150000e+01],
       [  9.16200000e-03,   3.17160000e+01],
       [  8.33600000e-03,   3.23720000e+01],
       [  7.51000000e-03,   3.31010000e+01],
       [  6.81800000e-03,   3.37820000e+01],
       [  6.12600000e-03,   3.45390000e+01],
       [  5.54800000e-03,   3.52440000e+01],
       [  4.97100000e-03,   3.60300000e+01],
       [  4.49200000e-03,   3.67590000e+01],
       [  4.01300000e-03,   3.75740000e+01],
       [  3.61900000e-03,   3.83280000e+01],
       [  3.22400000e-03,   3.91730000e+01],
       [  2.90000000e-03,   3.99510000e+01],
       [  2.57600000e-03,   4.08250000e+01],
       [  2.31200000e-03,   4.16270000e+01],
       [  2.04800000e-03,   4.25290000e+01],
       [  1.83400000e-03,   4.33550000e+01],
       [  1.61900000e-03,   4.42860000e+01],
       [  1.44600000e-03,   4.51340000e+01],
       [  1.27400000e-03,   4.60920000e+01],
       [  1.13500000e-03,   4.69620000e+01],
       [  9.96000000e-04,   4.79460000e+01],
       [  8.86000000e-04,   4.88350000e+01],
       [  7.75000000e-04,   4.98440000e+01],
       [  6.87000000e-04,   5.07540000e+01],
       [  5.99000000e-04,   5.17880000e+01],
       [  5.29000000e-04,   5.27160000e+01],
       [  4.60000000e-04,   5.37730000e+01],
       [  4.05000000e-04,   5.47170000e+01],
       [  3.50000000e-04,   5.57940000e+01],
       [  3.08000000e-04,   5.67520000e+01],
       [  2.65000000e-04,   5.78460000e+01],
       [  2.32000000e-04,   5.88160000e+01],
       [  1.99000000e-04,   5.99240000e+01],
       [  1.73000000e-04,   6.09020000e+01],
       [  1.48000000e-04,   6.20210000e+01],
       [  1.28000000e-04,   6.30040000e+01],
       [  1.08000000e-04,   6.41300000e+01],
       [  9.30000000e-05,   6.51150000e+01],
       [  7.80000000e-05,   6.62450000e+01],
       [  6.70000000e-05,   6.72430000e+01],
       [  5.50000000e-05,   6.83920000e+01],
       [  4.60000000e-05,   6.94400000e+01],
       [  3.70000000e-05,   7.06570000e+01],
       [  3.00000000e-05,   7.18120000e+01],
       [  2.20000000e-05,   7.31800000e+01],
       [  1.60000000e-05,   7.45940000e+01],
       [  1.00000000e-05,   7.63570000e+01],
       [  5.00000000e-06,   7.81460000e+01],
       [  0.00000000e+00,   8.05810000e+01]])
g5_eta_z = np.array([[  1.00000000e+00,  -6.00000000e-03],
       [  9.92500000e-01,   5.80000000e-02],
       [  9.84999000e-01,   1.23000000e-01],
       [  9.77456000e-01,   1.89000000e-01],
       [  9.69913000e-01,   2.54000000e-01],
       [  9.62370000e-01,   3.20000000e-01],
       [  9.54828000e-01,   3.87000000e-01],
       [  9.47285000e-01,   4.54000000e-01],
       [  9.39743000e-01,   5.21000000e-01],
       [  9.32200000e-01,   5.89000000e-01],
       [  9.24658000e-01,   6.57000000e-01],
       [  9.17116000e-01,   7.26000000e-01],
       [  9.09573000e-01,   7.95000000e-01],
       [  9.02031000e-01,   8.64000000e-01],
       [  8.94489000e-01,   9.34000000e-01],
       [  8.86948000e-01,   1.00400000e+00],
       [  8.79406000e-01,   1.07500000e+00],
       [  8.71864000e-01,   1.14600000e+00],
       [  8.64323000e-01,   1.21800000e+00],
       [  8.56781000e-01,   1.29000000e+00],
       [  8.49239000e-01,   1.36300000e+00],
       [  8.41698000e-01,   1.43600000e+00],
       [  8.34157000e-01,   1.51000000e+00],
       [  8.26616000e-01,   1.58400000e+00],
       [  8.19075000e-01,   1.65900000e+00],
       [  8.09021000e-01,   1.75900000e+00],
       [  7.98967000e-01,   1.86000000e+00],
       [  7.86400000e-01,   1.98800000e+00],
       [  7.73832000e-01,   2.11800000e+00],
       [  7.61265000e-01,   2.24900000e+00],
       [  7.48698000e-01,   2.38200000e+00],
       [  7.36134000e-01,   2.51700000e+00],
       [  7.23570000e-01,   2.65400000e+00],
       [  7.11006000e-01,   2.79200000e+00],
       [  6.98442000e-01,   2.93200000e+00],
       [  6.85878000e-01,   3.07400000e+00],
       [  6.73314000e-01,   3.21900000e+00],
       [  6.54471000e-01,   3.43900000e+00],
       [  6.35628000e-01,   3.66500000e+00],
       [  6.16790000e-01,   3.89600000e+00],
       [  5.97953000e-01,   4.13200000e+00],
       [  5.79115000e-01,   4.37500000e+00],
       [  5.60278000e-01,   4.62300000e+00],
       [  5.41449000e-01,   4.87900000e+00],
       [  5.22620000e-01,   5.14200000e+00],
       [  5.03795000e-01,   5.41300000e+00],
       [  4.84970000e-01,   5.69200000e+00],
       [  4.66153000e-01,   5.98000000e+00],
       [  4.47337000e-01,   6.27700000e+00],
       [  4.28528000e-01,   6.58500000e+00],
       [  4.09720000e-01,   6.90500000e+00],
       [  3.90927000e-01,   7.23700000e+00],
       [  3.72133000e-01,   7.58200000e+00],
       [  3.53349000e-01,   7.94300000e+00],
       [  3.34566000e-01,   8.32000000e+00],
       [  3.09854000e-01,   8.84600000e+00],
       [  2.85142000e-01,   9.40900000e+00],
       [  2.63587000e-01,   9.93600000e+00],
       [  2.42032000e-01,   1.05040000e+01],
       [  2.23772000e-01,   1.10210000e+01],
       [  2.05513000e-01,   1.15780000e+01],
       [  1.90061000e-01,   1.20860000e+01],
       [  1.74608000e-01,   1.26330000e+01],
       [  1.61513000e-01,   1.31340000e+01],
       [  1.48418000e-01,   1.36740000e+01],
       [  1.37287000e-01,   1.41700000e+01],
       [  1.26157000e-01,   1.47060000e+01],
       [  1.16695000e-01,   1.51980000e+01],
       [  1.07233000e-01,   1.57310000e+01],
       [  9.91910000e-02,   1.62220000e+01],
       [  9.11490000e-02,   1.67530000e+01],
       [  8.43130000e-02,   1.72430000e+01],
       [  7.74770000e-02,   1.77730000e+01],
       [  7.16000000e-02,   1.82690000e+01],
       [  6.57230000e-02,   1.88070000e+01],
       [  6.06820000e-02,   1.93090000e+01],
       [  5.56410000e-02,   1.98550000e+01],
       [  5.13260000e-02,   2.03640000e+01],
       [  4.70110000e-02,   2.09200000e+01],
       [  4.33260000e-02,   2.14380000e+01],
       [  3.96410000e-02,   2.20040000e+01],
       [  3.64990000e-02,   2.25310000e+01],
       [  3.33580000e-02,   2.31080000e+01],
       [  3.06730000e-02,   2.36480000e+01],
       [  2.79870000e-02,   2.42400000e+01],
       [  2.56990000e-02,   2.47940000e+01],
       [  2.34100000e-02,   2.54020000e+01],
       [  2.14670000e-02,   2.59710000e+01],
       [  1.95230000e-02,   2.65960000e+01],
       [  1.78780000e-02,   2.71800000e+01],
       [  1.62320000e-02,   2.78240000e+01],
       [  1.48440000e-02,   2.84230000e+01],
       [  1.34550000e-02,   2.90850000e+01],
       [  1.22870000e-02,   2.97010000e+01],
       [  1.11200000e-02,   3.03820000e+01],
       [  1.01410000e-02,   3.10150000e+01],
       [  9.16200000e-03,   3.17160000e+01],
       [  8.33600000e-03,   3.23720000e+01],
       [  7.51000000e-03,   3.31010000e+01],
       [  6.81800000e-03,   3.37820000e+01],
       [  6.12600000e-03,   3.45390000e+01],
       [  5.54800000e-03,   3.52440000e+01],
       [  4.97100000e-03,   3.60300000e+01],
       [  4.49200000e-03,   3.67590000e+01],
       [  4.01300000e-03,   3.75740000e+01],
       [  3.61900000e-03,   3.83280000e+01],
       [  3.22400000e-03,   3.91730000e+01],
       [  2.90000000e-03,   3.99510000e+01],
       [  2.57600000e-03,   4.08250000e+01],
       [  2.31200000e-03,   4.16270000e+01],
       [  2.04800000e-03,   4.25290000e+01],
       [  1.83400000e-03,   4.33550000e+01],
       [  1.61900000e-03,   4.42860000e+01],
       [  1.44600000e-03,   4.51340000e+01],
       [  1.27400000e-03,   4.60920000e+01],
       [  1.13500000e-03,   4.69620000e+01],
       [  9.96000000e-04,   4.79460000e+01],
       [  8.86000000e-04,   4.88350000e+01],
       [  7.75000000e-04,   4.98440000e+01],
       [  6.87000000e-04,   5.07540000e+01],
       [  5.99000000e-04,   5.17880000e+01],
       [  5.29000000e-04,   5.27160000e+01],
       [  4.60000000e-04,   5.37730000e+01],
       [  4.05000000e-04,   5.47170000e+01],
       [  3.50000000e-04,   5.57940000e+01],
       [  3.08000000e-04,   5.67520000e+01],
       [  2.65000000e-04,   5.78460000e+01],
       [  2.32000000e-04,   5.88160000e+01],
       [  1.99000000e-04,   5.99240000e+01],
       [  1.73000000e-04,   6.09020000e+01],
       [  1.48000000e-04,   6.20210000e+01],
       [  1.28000000e-04,   6.30040000e+01],
       [  1.08000000e-04,   6.41300000e+01],
       [  9.30000000e-05,   6.51150000e+01],
       [  7.80000000e-05,   6.62450000e+01],
       [  6.70000000e-05,   6.72430000e+01],
       [  5.50000000e-05,   6.83920000e+01],
       [  4.60000000e-05,   6.94400000e+01],
       [  3.70000000e-05,   7.06570000e+01],
       [  3.00000000e-05,   7.18120000e+01],
       [  2.20000000e-05,   7.31800000e+01],
       [  1.60000000e-05,   7.45940000e+01],
       [  1.00000000e-05,   7.63570000e+01],
       [  5.00000000e-06,   7.81460000e+01],
       [  0.00000000e+00,   8.05810000e+01]])

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
    def __init__(self, path, mode = 'r'):
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
                       eta_z = g5_eta_z
                       nlay = 72
                       nlay_stag = 73
                   elif name[:2] == 'G4':
                       eta_z = g5_eta_z
                       nlay = 55
                       nlay_stag = 56
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
        self.createVariable('layer', 'f', ('layer',), values = eta_z[1::2, 0], units = 'hybrid pressure-sigma', bounds = 'layer_bounds')
        self.createVariable('layer_bounds', 'f', ('layer','nv'), values = eta_z[::2, 0].repeat(2, 0)[1:-1].reshape(-1, 2), units = 'hybrid pressure-sigma')
        self.createVariable('layer_height', 'f', ('layer',), values = eta_z[1::2, 1], units = 'km', bounds = 'layer_height_bounds')
        self.createVariable('layer_height_bounds', 'f', ('layer','nv'), values = eta_z[::2, 1].repeat(2, 0)[1:-1].reshape(-1, 2), units = 'km')
        
        self.createVariable('time', 'f', ('time',), values = np.array([(datetime(y, m, d, H, M, S) - datetime(1985, 1, 1, 0, 0, 0)).total_seconds() / 3600. for y, m, d, H, M, S in zip(years, months, days, hours, minutes, seconds)]), units = 'hours since 1985-01-01 00:00:00 UTC', base_units = 'hours since 1985-01-01 00:00:00 UTC', long_name = 'time', standard_name = 'time')
        self.createVariable('longitude', 'f', ('longitude',), values = longitude, units = 'degrees_east', standard_name = 'longitude', bounds = 'longitude_bounds')
        self.createVariable('longitude_bounds', 'f', ('longitude','nv'), values = longitude_bounds, units = 'degrees_east', standard_name = 'longitude')
        self.createVariable('latitude', 'f', ('latitude',), values = latitude, units = 'degrees_north', standard_name = 'latitude', bounds = 'latitude_bounds')
        self.createVariable('latitude_bounds', 'f', ('latitude','nv'), values = latitude_bounds, units = 'degrees_north', standard_name = 'latitude')

if __name__ == '__main__':
    path = sys.argv[1]
    f = geos(path)
