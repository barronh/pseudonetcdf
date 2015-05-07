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
            
            self._data = np.concatenate(pieces, axis = 0).view(datat)['data'] * self.scale
        return self._data.__getitem__(k)
            
    def ncattrs(self):
        return tuple(self.__dict__.keys())
    
class bpch2(PseudoNetCDFFile):
    def __init__(self, path):
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
            tmpvariables[key] = gcvar(key, self)
        longitudes = []
        latitudes = []
        levels = []
        for key in outpos:
            tmpvar = tmpvariables[key]
            longitudes.append(tmpvar.nlongitudes)
            latitudes.append(tmpvar.nlatitudes)
            levels.append(tmpvar.nlevels)
            self.createVariable(tmpvar.category + '_'+ tmpvar.shortname, tmpvar.dtype.char, ('time', 'layer%d' % tmpvar.nlevels, 'latitude', 'longitude'), values = tmpvar)
        self.createDimension('time', max([len(pos) for pos in outpos.values()]))
        self.createDimension('layer', max(levels))
        for layer in levels:
            self.createDimension('layer%d' % layer, layer)
        self.createDimension('latitude', max(latitudes)) 
        self.createDimension('longitude', max(longitudes)) 
        
    def _gettracerinfo(self, path):
        tpath = os.path.join(os.path.dirname(path) + 'tracerinfo.dat')
        self._tdata = np.recfromtxt(tpath, dtype=None, comments = '#', names=['shortname','fullname','molwt', 'ncarbon', 'tracer', 'scale', 'unit'], delimiter=[9,30,10,3,9,10,41], autostrip = True) 

    def _getdiaginfo(self, path):
        dpath = os.path.join(os.path.dirname(path) + 'diaginfo.dat')
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