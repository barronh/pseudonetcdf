import os
from collections import OrderedDict
import numpy as np
from PseudoNetCDF.sci_var import PseudoNetCDFFile

_datablock_header_type = np.dtype([('bpad1', '>i4'),
                                    ('modelname', 'S20'),
                                    ('modelres', '>f4', 2),
                                    ('halfpolar', '>i4'),
                                    ('center180', '>i4'),
                                    ('epad1', '>i4'),
                                    ('bpad2', '>i4'),
                                    ('category', 'S40'),
                                    ('tracer', '>i4'),
                                    ('unit', 'S40'),
                                    ('tau0', '>f8'),
                                    ('tau1', '>f8'),
                                    ('reserved', 'S40'),
                                    ('dim', '>i4', 3),
                                    ('start', '>i4', 3),
                                    ('skip', '>i4'),
                                    ('epad2', '>i4')])
_hdr_size = _datablock_header_type.itemsize

def add_vert(self, key):
    from _vertcoord import geos_etai_pressure, geos_etam_pressure, geos_hyam, geos_hyai, geos_hybm, geos_hybi
    if key == 'hyai':
        data = self.Ap
        dims = ('layer_bounds', )
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(units = "hPa", long_name = "hybrid A coefficient at layer interfaces", note = "unit consistent with GEOS-Chem pressure outputs")
    elif key == 'hyam':
        data = geos_hyam[self.vertgrid]
        dims = ('layer', )
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(units = "hPa", long_name = "hybrid B coefficient at layer midpoints", note = "unit consistent with GEOS-Chem pressure outputs")
    elif key == 'hybi':
        data = self.Bp
        dims = ('layer_bounds',)
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(units = "1", long_name = "hybrid B coefficient at layer interfaces")
    elif key == 'hybm':
        data = geos_hybm[self.vertgrid]
        dims = ('layer', )
        dtype = data.dtype.char
        if dims[0] not in self.dimensions:
            self.createDimension(dims[0], data.size)
        kwds = dict(units = "1", long_name = "hybrid B coefficient at layer midpoints")
    tmpvar = self.createVariable(key, 'f', dims)
    
    for p, v in kwds.items():
        setattr(tmpvar, p, v)
    
    tmpvar[:] = data

def add_lat(self, key, example):
    yres = self.modelres[1]
    if self.halfpolar == 1:
        data = np.concatenate([[-90.], np.arange(-90. + yres / 2., 90., yres), [90.]])
    else:
        data = np.arange(-90, 90 + yres, yres)
    dims = ('latitude',)
    dtype = 'f'
    kwds = dict(standard_name = "latitude", long_name = "latitude", units = "degrees_north", base_units = "degrees_north", axis = "Y")
    if key == 'latitude':
        data = data[:-1] + np.diff(data) / 2.
        kwds['bounds'] = 'latitude_bounds'
    else:
        dims += ('nv',)
        data = data.repeat(2,0)[1:-1].reshape(-1, 2)
    sj = getattr(example, 'STARTJ', 0)
    data = data[sj:sj + example.shape[2]]
    var = self.createVariable(key, 'f', dims)
    var[:] = data
    for p, v in kwds.items():
        setattr(self, p, v)

def add_lon(self, key, example):
    xres = self.modelres[0]
    i = np.arange(0, 360 + xres, xres)
    data = i - (180 + xres / 2. * self.center180)
    dims = ('longitude',)
    dtype = 'f'
    kwds = dict(standard_name = "longitude", long_name = "longitude", units = "degrees_east", base_units = "degrees_east", axis = "X")
    if key == 'longitude':
        data = data[:-1] + np.diff(data) / 2.
        kwds['bounds'] = 'longitude_bounds'
    else:
        dims += ('nv',)
        data = data.repeat(2,0)[1:-1].reshape(-1, 2)
    si = getattr(example, 'STARTI', 0)
    data = data[si:si + example.shape[3]]
    var = self.createVariable(key, 'f', dims)
    var[:] = data
    for p, v in kwds.items():
        setattr(self, p, v)

class gcvar(object):
    def __init__(self, key, parent):
        self._key = key
        self._parent = parent
        self._data = None
        start, end, dim = self._parent._outpos[self._key].values()[0]
        self.nlevels = dim[2]
        self.nlatitudes = dim[1]
        self.nlongitudes = dim[0]
        #self.dimensions = ('time', 'layer%d' % self.nlevels, 'latitude', 'longitude')
        self._header = self._parent._data[start:start+_hdr_size].view(_datablock_header_type)
        self.category = self._header['category'][0].strip()
        self.tracer = self._header['tracer'][0]
        self.catoffset = [row['offset'] for row in self._parent._ddata if row['category'] == self.category][0]
        self.cattracer = self.catoffset + self.tracer
        props = ([row for row in self._parent._tdata if row['tracer'] == self.cattracer] + [row for row in self._parent._tdata if row['tracer'] == self.tracer])[0]
        for pk in props.dtype.names:
            setattr(self, pk, props[pk])
        setattr(self, 'units', getattr(self, 'unit'))
    def __getattr__(self, k):
        try:
            return object.__getattribute__(self, k)
        except:
            return object.__getattribute__(self[:], k)
    
    def __getitem__(self, k):
        if self._data is None:
            pieces = []
            for (tstart, tend), (start, end, dim) in self._parent._outpos[self._key].iteritems():
                nelem = end - start - _hdr_size
                datat = np.dtype([('header', _datablock_header_type), ('bpad', '>i'), ('data', '>f', dim[::-1]), ('epad', '>i')])
                
                pieces.append(self._parent._data[start:end])
            tmpdata = np.concatenate(pieces, axis = 0).view(datat)
            self._data = tmpdata['data'] * self.scale
            self._tau0 = tmpdata['header']['tau0']
            self._tau1 = tmpdata['header']['tau1']
        return self._data.__getitem__(k)
            
    def ncattrs(self):
        return tuple(self.__dict__.keys())
    
class bpch2(PseudoNetCDFFile):
    def __init__(self, path, nogroup = False, vertgrid = 'GEOS-5-REDUCED'):
        from _vertcoord import geos_hyai, geos_hybi
        self.vertgrid = vertgrid
        # Ap [hPa]
        self.Ap = geos_hyai[self.vertgrid]

        # Bp [unitless]
        self.Bp = geos_hybi[self.vertgrid]

        self._gettracerinfo(path)
        self._getdiaginfo(path)
        self._data = data = np.memmap(path, mode = 'r', dtype = 'uint8')
        header = data[:136]
        datastart = offset = 136

        outpos = self._outpos = OrderedDict()
        while offset < data.size:
            tmp_hdr = data[offset:offset+_hdr_size].view(_datablock_header_type)
            key = tmp_hdr['category'][0].strip() + '_' + str(tmp_hdr['tracer'][0])
            tau0, tau1 = tmp_hdr['tau0'][0], tmp_hdr['tau1'][0]
            dim = tuple(tmp_hdr['dim'][0].tolist())
            skip = tmp_hdr['skip'][0]
            thispos = outpos.setdefault(key, OrderedDict())[(tau0, tau1)] = offset, offset + _hdr_size + skip, dim
            offset += skip + _hdr_size
                
        tmpvariables = {}
        for key in outpos:
            tmpvar = tmpvariables[key] = gcvar(key, self)
        longitudes = []
        latitudes = []
        levels = []
        self.modelres = tmpvar._header['modelres'][0]
        self.halfpolar = tmpvar._header['halfpolar'][0]
        self.center180 = tmpvar._header['center180'][0]
        for key in outpos:
            tmpvar = tmpvariables[key]
            longitudes.append(tmpvar.nlongitudes)
            latitudes.append(tmpvar.nlatitudes)
            levels.append(tmpvar.nlevels)
            if nogroup:
                tmpkey = tmpvar.shortname
            else:
                tmpkey = tmpvar.category + '_'+ tmpvar.shortname
            var = self.createVariable(tmpkey, tmpvar.dtype.char, ('time', 'layer%d' % tmpvar.nlevels, 'latitude', 'longitude'), values = tmpvar)
            for k in tmpvar.ncattrs():
                setattr(var, k, getattr(tmpvar, k))
        tmpvar = tmpvariables.values()[0]
        self.createDimension('time', max([len(pos) for pos in outpos.values()]))
        self.createDimension('layer', max(levels))
        for layer in levels:
            self.createDimension('layer%d' % layer, layer)
        self.createDimension('latitude', max(latitudes)) 
        self.createDimension('longitude', max(longitudes)) 
        self.createDimension('nv', 2)
        add_lat(self, 'latitude', tmpvar)
        add_lat(self, 'latitude_bounds', tmpvar)
        add_lon(self, 'longitude', tmpvar)
        add_lon(self, 'longitude_bounds', tmpvar)
        add_vert(self, 'hyam')
        add_vert(self, 'hyai')
        add_vert(self, 'hybm')
        add_vert(self, 'hybi')
        for k, v in [('tau0', tmpvar._tau0), ('time', tmpvar._tau0), ('tau1', tmpvar._tau1)]:
            tvar = self.createVariable(k, 'i', ('time',))
            tvar.units = 'hours since 1985-01-01 00:00:00 UTC'
            tvar[:] = v[:]
        
    def _gettracerinfo(self, path):
        tpath = os.path.join(os.path.dirname(path), 'tracerinfo.dat')
        if not os.path.exists(tpath):
            tpath = 'tracerinfo.dat'
        self._tdata = np.recfromtxt(tpath, dtype=None, comments = '#', names=['shortname','fullname','kgpermole', 'carbon', 'tracer', 'scale', 'unit'], delimiter=[9,30,10,3,9,10,41], autostrip = True) 

    def _getdiaginfo(self, path):
        dpath = os.path.join(os.path.dirname(path), 'diaginfo.dat')
        if not os.path.exists(dpath):
            dpath = 'diaginfo.dat'
        self._ddata = np.recfromtxt(dpath, dtype=None, comments = '#', names=['offset','category','comment'], delimiter=[9,40,100], autostrip = True)
        
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

if __name__ == '__main__':
    f = bpch2('/Users/barronh/Development/pseudonetcdf/src/PseudoNetCDF/testcase/geoschemfiles/test.bpch')
    var = f.variables['IJ-AVG-$_Ox']
    print var[:].mean()