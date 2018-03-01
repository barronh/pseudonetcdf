from __future__ import print_function
import os
from collections import OrderedDict
import numpy as np
from PseudoNetCDF.sci_var import PseudoNetCDFVariables
from ._bpch import _general_header_type, bpch as oldbpch, bpch_base
import unittest


_datablock_header_type = np.dtype([('bpad1', '>i4'),
                                   ('modelname', 'S20'),
                                   ('modelres', '>f4', 2),
                                   ('halfpolar', '>i4'),
                                   ('center180', '>i4'),
                                   ('epad1', '>i4'),
                                   ('bpad2', '>i4'),
                                   ('category', 'S40'),
                                   ('tracerid', '>i4'),
                                   ('base_units', 'S40'),
                                   ('tau0', '>f8'),
                                   ('tau1', '>f8'),
                                   ('reserved', 'S40'),
                                   ('dim', '>i4', 3),
                                   ('start', '>i4', 3),
                                   ('skip', '>i4'),
                                   ('epad2', '>i4')])
_hdr_size = _datablock_header_type.itemsize


def add_vert(self, key):
    from ._vertcoord import geos_etai_pressure, geos_etam_pressure
    from ._vertcoord import geos_hyam, geos_hybm
    if key == 'hyai':
        data = self.Ap
        dims = ('layer_bounds', )
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(units="hPa",
                    long_name="hybrid A coefficient at layer interfaces",
                    note="unit consistent with GEOS-Chem pressure outputs")
    elif key == 'hyam':
        data = geos_hyam[self.vertgrid]
        dims = ('layer',)
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(units="hPa",
                    long_name="hybrid A coefficient at layer midpoints",
                    note="unit consistent with GEOS-Chem pressure outputs")
    elif key == 'hybi':
        data = self.Bp
        dims = ('layer_bounds',)
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(
            units="1", long_name="hybrid B coefficient at layer interfaces")
    elif key == 'etai_pressure':
        data = geos_etai_pressure[self.vertgrid]
        dtype = data.dtype.char
        dims = ('layer_bounds',)
        kwds = dict(units='hPa', long_name='eta levels at interfaces')
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
    elif key == 'etam_pressure':
        data = geos_etam_pressure[self.vertgrid]
        dtype = data.dtype.char
        dims = ('layer',)
        kwds = dict(units='hPa', long_name='eta levels at mid points')
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
    elif key == 'hybm':
        data = geos_hybm[self.vertgrid]
        dims = ('layer', )
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(
            units="1", long_name="hybrid B coefficient at layer midpoints")
    tmpvar = self.createVariable(key, dtype, dims)
    nlay = tmpvar.size

    for p, v in kwds.items():
        setattr(tmpvar, p, v)

    tmpvar[:] = data[0:nlay]


def add_lat(self, key, example):
    yres = self.modelres[1]
    if self.halfpolar == 1:
        data = np.concatenate(
            [[-90.], np.arange(-90. + yres / 2., 90., yres), [90.]])
    else:
        data = np.arange(-90, 90 + yres, yres)
    dims = ('latitude',)
    dtype = 'f'
    kwds = dict(standard_name="latitude", long_name="latitude",
                units="degrees_north", base_units="degrees_north", axis="Y")
    if key == 'latitude':
        data = data[:-1] + np.diff(data) / 2.
        kwds['bounds'] = 'latitude_bounds'
    else:
        dims += ('nv',)
        data = data.repeat(2, 0)[1:-1].reshape(-1, 2)
    sj = getattr(example, 'STARTJ', 0)
    data = data[sj:sj + example.shape[2]]
    var = self.createVariable(key, dtype, dims)
    var[:] = data
    for p, v in kwds.items():
        setattr(var, p, v)


def add_lon(self, key, example):
    xres = self.modelres[0]
    i = np.arange(0, 360 + xres, xres)
    data = i - (180 + xres / 2. * self.center180)
    dims = ('longitude',)
    dtype = 'f'
    kwds = dict(standard_name="longitude", long_name="longitude",
                units="degrees_east", base_units="degrees_east", axis="X")
    if key == 'longitude':
        data = data[:-1] + np.diff(data) / 2.
        kwds['bounds'] = 'longitude_bounds'
    else:
        dims += ('nv',)
        data = data.repeat(2, 0)[1:-1].reshape(-1, 2)
    si = getattr(example, 'STARTI', 0)
    data = data[si:si + example.shape[3]]
    var = self.createVariable(key, dtype, dims)
    var[:] = data
    for p, v in kwds.items():
        setattr(var, p, v)


class gcvar(object):
    def __init__(self, key, parent):
        self._key = key
        self._name = key
        self._parent = parent
        self._data = None
        start, end, dim = list(self._parent._outpos[self._key].values())[0]
        ntimes = len(self._parent._outpos[self._key])
        nlevels = dim[2]
        nlatitudes = dim[1]
        nlongitudes = dim[0]
        timedim = 'time'
        laydim = 'layer%d' % nlevels
        self.dimensions = (timedim, laydim, 'latitude', 'longitude')
        self.shape = (ntimes, nlevels, nlatitudes, nlongitudes)
        dht = _datablock_header_type
        self._header = self._parent._data[start:start +
                                          _hdr_size].view(dht)
        self.category = self._header['category'][0].strip()
        self.tracerid = self._header['tracerid'][0]
        self.base_units = self._header['base_units'][0]
        self.catoffset = [row['offset']
                          for row in self._parent._ddata
                          if row['category'] == self.category][0]
        self.noscale = self._parent.noscale
        STARTK, STARTJ, STARTI = self._header['start'][0][::-1] - 1
        self.STARTK, self.STARTJ, self.STARTI = STARTK, STARTJ, STARTI
        if hasattr(self.category, 'decode'):
            self.category = self.category.decode()
        self.cattracerid = self.catoffset + self.tracerid
        props = ([row for row in self._parent._tdata
                  if row['tracerid'] == self.cattracerid] +
                 [row for row in self._parent._tdata
                  if row['tracerid'] == self.tracerid])[0]
        for pk in props.dtype.names:
            if pk == 'tracerid':
                continue
            pv = props[pk]
            if hasattr(pv, 'decode'):
                pv = pv.decode()
            setattr(self, pk, pv)

    def __getattr__(self, k):
        try:
            return object.__getattribute__(self, k)
        except Exception:
            return object.__getattribute__(self[:], k)

    def __getitem__(self, k):
        if self._data is None:
            pieces = []
            startenddims = self._parent._outpos[self._key]
            for (tstart, tend), (start, end, dim) in startenddims.items():
                # nelem = end - start - _hdr_size
                datat = np.dtype([('header', _datablock_header_type),
                                  ('bpad', '>i'),
                                  ('data', '>f', dim[::-1]),
                                  ('epad', '>i')])

                pieces.append(self._parent._data[start:end])
            tmpdata = np.concatenate(pieces, axis=0).view(datat)
            if self.noscale:
                self._data = tmpdata['data']
            else:
                self._data = tmpdata['data'] * self.scale
            self._tau0 = tmpdata['header']['tau0']
            self._tau1 = tmpdata['header']['tau1']
        return self._data.__getitem__(k)

    def ncattrs(self):
        return tuple(self.__dict__.keys())


class bpch2(bpch_base):
    @classmethod
    def isMine(cls, path):
        isbpch = oldbpch.isMine(path)
        if not isbpch:
            return False

        try:
            oldbpch(path)
            return False
        except Exception:
            return True

    def ij2ll(self, i, j):
        """
        i, j to center lon, lat
        """
        lat = np.asarray(self.variables['latitude'])
        lon = np.asarray(self.variables['longitude'])
        return lon[i], lat[j]

    def ll2ij(self, lon, lat):
        """
        lon and lat may be scalars or vectors, but matrices will crash
        """
        late = np.asarray(self.variables['latitude_bounds'])
        lone = np.asarray(self.variables['longitude_bounds'])
        inlon = (lon >= lone[:, 0, None]) & (lon < lone[:, 1, None])
        inlat = (lat >= late[:, 0, None]) & (lat < late[:, 1, None])
        i = inlon.argmax(0)
        j = inlat.argmax(0)
        return i, j

    def __init__(self, path, nogroup=False, noscale=False,
                 vertgrid='GEOS-5-REDUCED'):
        self.noscale = noscale
        self._nogroup = nogroup
        from ._vertcoord import geos_hyai, geos_hybi
        self.vertgrid = vertgrid
        # Ap [hPa]
        self.Ap = geos_hyai[self.vertgrid]

        # Bp [unitless]
        self.Bp = geos_hybi[self.vertgrid]

        self._gettracerinfo(path)
        self._getdiaginfo(path)
        self._data = data = np.memmap(path, mode='r', dtype='uint8')
        header = data[:136].copy().view(_general_header_type)
        self.ftype = header[0][1]
        self.toptitle = header[0][4]

        offset = 136

        outpos = self._outpos = OrderedDict()

        while offset < data.size:
            tmp_hdr = data[offset:offset +
                           _hdr_size].view(_datablock_header_type)
            cat = tmp_hdr['category'][0].strip()
            if hasattr(cat, 'decode'):
                cat = cat.decode()
            trac = str(tmp_hdr['tracerid'][0])
            if hasattr(trac, 'decode'):
                trac = trac.decode()

            key = cat + '_' + trac
            tau0, tau1 = tmp_hdr['tau0'][0], tmp_hdr['tau1'][0]
            dim = tuple(tmp_hdr['dim'][0].tolist())
            skip = tmp_hdr['skip'][0]
            outpos.setdefault(key, OrderedDict())[(
                tau0, tau1)] = offset, offset + _hdr_size + skip, dim
            offset += skip + _hdr_size

        modelname, modelres, halfpolar, center180 = tmp_hdr[0].tolist()[1:5]
        self.modelname = modelname
        self.modelres = modelres
        self.halfpolar = halfpolar
        self.center180 = center180

        tmpvariables = self._gcvars = OrderedDict()
        for key in outpos:
            tmpvar = gcvar(key, self)
            if nogroup is True:
                tmpkey = str(tmpvar.shortname)
            elif nogroup is False:
                tmpkey = tmpvar.category + str('_') + tmpvar.shortname
            elif tmpvar.category in nogroup:
                tmpkey = str(tmpvar.shortname)
            else:
                tmpkey = tmpvar.category + str('_') + tmpvar.shortname

            tmpvariables[tmpkey] = tmpvar
        self.modelres = tmpvar._header['modelres'][0]
        self.halfpolar = tmpvar._header['halfpolar'][0]
        self.center180 = tmpvar._header['center180'][0]
        ntimes = [tmpvar.shape[0] for tmpkey, tmpvar in tmpvariables.items()]
        levels = [tmpvar.shape[1] for tmpkey, tmpvar in tmpvariables.items()]
        latitudes = [tmpvar.shape[2]
                     for tmpkey, tmpvar in tmpvariables.items()]
        longitudes = [tmpvar.shape[3]
                      for tmpkey, tmpvar in tmpvariables.items()]
        self._maxntimes = maxntimes = max(ntimes)
        self.createDimension('time', maxntimes)
        for time in set(ntimes):
            if time != maxntimes:
                self.createDimension('time%d' % time, time)
        self.createDimension('layer', min(max(levels), self.Ap.size - 1))
        for layer in set(levels):
            self.createDimension('layer%d' % layer, layer)
        self.createDimension('latitude', max(latitudes))
        self.createDimension('longitude', max(longitudes))
        self.createDimension('nv', 2)
        self.variables = PseudoNetCDFVariables(
            self._addvar, list(tmpvariables))
        key = list(tmpvariables)[0]
        tmpvar = self.variables[key]
        add_lat(self, 'latitude', tmpvar)
        add_lat(self, 'latitude_bounds', tmpvar)
        add_lon(self, 'longitude', tmpvar)
        add_lon(self, 'longitude_bounds', tmpvar)
        add_vert(self, 'hyam')
        add_vert(self, 'hyai')
        add_vert(self, 'hybm')
        add_vert(self, 'hybi')
        add_vert(self, 'etai_pressure')
        add_vert(self, 'etam_pressure')
        for k, v in [('tau0', tmpvar._tau0),
                     ('time', tmpvar._tau0),
                     ('tau1', tmpvar._tau1)]:
            tvar = self.createVariable(k, 'i', ('time',))
            tvar.units = 'hours since 1985-01-01 00:00:00 UTC'
            tvar[:] = v[:]

    def _addvar(self, key):
        tmpvar = self._gcvars[key]
        if tmpvar.shape[0] == self._maxntimes:
            timedim = 'time'
        else:
            timedim = 'time%d' % tmpvar.shape[0]
        laydim = 'layer%d' % tmpvar.shape[1]
        outdims = (timedim, laydim, 'latitude', 'longitude')
        var = self.createVariable(
            key, tmpvar.dtype.char, outdims, values=tmpvar)
        for k in tmpvar.ncattrs():
            if k in ('_name', 'shape', 'dimensions'):
                continue
            setattr(var, k, getattr(tmpvar, k))
        return var

    def _gettracerinfo(self, path):
        tpath = os.path.join(os.path.dirname(path), 'tracerinfo.dat')
        if not os.path.exists(tpath):
            tpath = 'tracerinfo.dat'
        self._tdata = np.recfromtxt(tpath, dtype=None, comments='#', names=[
                                    'shortname', 'fullname', 'kgpermole',
                                    'carbon', 'tracerid', 'scale', 'units'],
                                    delimiter=[9, 30, 10, 3, 9, 10, 41],
                                    autostrip=True)

    def _getdiaginfo(self, path):
        dpath = os.path.join(os.path.dirname(path), 'diaginfo.dat')
        if not os.path.exists(dpath):
            dpath = 'diaginfo.dat'
        self._ddata = np.recfromtxt(dpath, dtype=None, comments='#', names=[
                                    'offset', 'category', 'comment'],
                                    delimiter=[9, 40, 100], autostrip=True)

# OFFSET    (I8 )  Constant to add to tracer numbers in order to distinguish
#                  for the given diagnostic category, as stored in file
#                  "tracerinfo.dat".  OFFSET may be up to 8 digits long.
#  --       (1X )  1-character spacer
# CATEGORY  (A40)  Category name for CTM diagnostics. NOTE: The category name
#                  can be up to 40 chars long, but historically the GEOS-CHEM
#                  and GISS models have used an 8-character category name.
# COMMENT   (A  )  Descriptive comment string
#
#  --       (1X )  1-character spacer


class TestMemmaps(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF.testcase import geoschemfiles_paths
        self.bpchpath = geoschemfiles_paths['bpch']

    def testNCF2BPCH(self):
        bpchfile = bpch2(self.bpchpath, noscale=True)
        from PseudoNetCDF.pncgen import pncgen
        pncgen(bpchfile, self.bpchpath + '.check', inmode='r',
               outmode='w', format='bpch', verbose=0)
        orig = open(self.bpchpath, 'rb').read()
        new = open(self.bpchpath + '.check', 'rb').read()
        assert(orig == new)
        os.remove(self.bpchpath+'.check')

    def testBPCH2(self):
        bpchfile = bpch2(self.bpchpath)
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        ALD2_check = np.array([1.60520077e-02, 1.82803553e-02, 2.00258084e-02,
                               2.01461259e-02, 1.84865110e-02, 2.49667447e-02,
                               2.73083989e-02, 2.87465211e-02, 2.89694592e-02,
                               2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
                               3.22086290e-02, 3.35262120e-02, 3.41329686e-02,
                               3.05218045e-02, 3.30278911e-02, 3.58164124e-02,
                               3.93186994e-02, 4.15412188e-02, 1.60520077e-02,
                               1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
                               1.84865110e-02, 2.49667447e-02, 2.73083989e-02,
                               2.87465211e-02, 2.89694592e-02, 2.87686456e-02,
                               2.87277419e-02, 3.08121163e-02, 3.22086290e-02,
                               3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
                               3.30278911e-02, 3.58164124e-02, 3.93186994e-02,
                               4.15412188e-02, 1.60520077e-02, 1.82803553e-02,
                               2.00258084e-02, 2.01461259e-02, 1.84865110e-02,
                               2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
                               2.89694592e-02, 2.87686456e-02, 2.87277419e-02,
                               3.08121163e-02, 3.22086290e-02, 3.35262120e-02,
                               3.41329686e-02, 3.05218045e-02, 3.30278911e-02,
                               3.58164124e-02, 3.93186994e-02,
                               4.15412188e-02]).reshape(ALD2.shape)
        np.testing.assert_allclose(ALD2, ALD2_check)
        hyai = np.array([0.0, 0.04804826, 6.593752,
                         13.1348, 19.61311, 26.09201, 32.57081, 38.98201,
                         45.33901, 51.69611, 58.05321, 64.36264, 70.62198,
                         78.83422, 89.09992, 99.36521, 109.1817, 118.9586,
                         128.6959, 142.91, 156.26, 169.609, 181.619, 193.097,
                         203.259, 212.15, 218.776, 223.898, 224.363, 216.865,
                         201.192, 176.93, 150.393, 127.837, 108.663, 92.36572,
                         78.51231, 56.38791, 40.17541, 28.36781, 19.7916,
                         9.292942, 4.076571, 1.65079, 0.6167791, 0.211349,
                         0.06600001, 0.01], 'f')
        np.testing.assert_allclose(bpchfile.variables['hyai'], hyai)
        hybi = np.array([1.0, 0.984952, 0.963406, 0.941865,
                         0.920387, 0.898908, 0.877429, 0.856018, 0.8346609,
                         0.8133039, 0.7919469, 0.7706375, 0.7493782, 0.721166,
                         0.6858999, 0.6506349, 0.6158184, 0.5810415, 0.5463042,
                         0.4945902, 0.4437402, 0.3928911, 0.3433811,
                         0.2944031, 0.2467411, 0.2003501, 0.1562241,
                         0.1136021, 0.06372006, 0.02801004, 0.006960025,
                         8.175413e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'f')
        np.testing.assert_allclose(bpchfile.variables['hybi'], hybi)
        etaip = np.array([1.01325000e+03, 9.98050662e+02,
                          9.82764882e+02, 9.67479511e+02, 9.52195238e+02,
                          9.36910541e+02, 9.21625744e+02, 9.06342248e+02,
                          8.91059167e+02, 8.75776287e+02, 8.60493406e+02,
                          8.45211087e+02, 8.29929441e+02, 8.09555669e+02,
                          7.84087994e+02, 7.58621022e+02, 7.33159694e+02,
                          7.07698900e+02, 6.82238631e+02, 6.44053520e+02,
                          6.05879758e+02, 5.67705907e+02, 5.29549900e+02,
                          4.91400941e+02, 4.53269420e+02, 4.15154739e+02,
                          3.77070069e+02, 3.39005328e+02, 2.88927351e+02,
                          2.45246173e+02, 2.08244245e+02, 1.76930008e+02,
                          1.50393000e+02, 1.27837000e+02, 1.08663000e+02,
                          9.23657200e+01, 7.85123100e+01, 5.63879100e+01,
                          4.01754100e+01, 2.83678100e+01, 1.97916000e+01,
                          9.29294200e+00, 4.07657100e+00, 1.65079000e+00,
                          6.16779100e-01, 2.11349000e-01, 6.60000100e-02,
                          1.00000000e-02])
        np.testing.assert_allclose(bpchfile.variables['etai_pressure'], etaip)

    def runTest(self):
        pass


if __name__ == '__main__':
    bpchpath = ('/Users/barronh/Development/pseudonetcdf/src/PseudoNetCDF/' +
                'testcase/geoschemfiles/test.bpch')
    f = bpch2(bpchpath)
    var = f.variables['IJ-AVG-$_Ox']
    print(var[:].mean())
