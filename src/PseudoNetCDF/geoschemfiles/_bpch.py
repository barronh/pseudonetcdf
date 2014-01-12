__all__ = ['bpch']
from matplotlib import use

import os
import gc
import re

# part of the default Python distribution
from collections import defaultdict
from warnings import warn
from StringIO import StringIO

# numpy is a very common installed library
import numpy as np
from numpy import ndarray, fromfile, memmap, dtype, arange, zeros, ceil, diff, concatenate, append, pi, sin

try:
    from PseudoNetCDF import PseudoNetCDFDimension, PseudoNetCDFVariable, PseudoNetCDFFile
except:
    warn('Using static PseudoNetCDF')
    # PseudoNetCDF is my own home grown
    # https://dawes.sph.unc.edu/trac/PseudoNetCDF
    #from PseudoNetCDF import PseudoNetCDFVariable, PseudoNetCDFFile
    class PseudoNetCDFDimension(object):
        """
        Dimension object responds like that of netcdf4-python
        """
        def __init__(self, group, name, size):
            self._len = int(size)
            self._unlimited = False
        def isunlimited(self):
            return self._unlimited
        def __len__(self):
            return self._len
        def setunlimited(self, unlimited):
            self._unlimited = unlimited

    class PseudoNetCDFFile(object):
        """
        PseudoNetCDFFile provides an interface and standard set of
        methods that a file should present to act like a netCDF file
        using the Scientific.IO.NetCDF.NetCDFFile interface.
        """
        def __new__(cls, *args, **kwds):
            new = object.__new__(cls, *args, **kwds)
            new.variables={}
            new.dimensions={}
            new._ncattrs = ()
            return new
    
        def __init__(self, *args, **properties):
            for k, v in properties.iteritems():
                setattr(self, k, v)

        def __setattr__(self, k, v):
            if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
                self._ncattrs += (k,)
            object.__setattr__(self, k, v)
        def createDimension(self,name,length):
            """
            name - string name for dimension
            length - maximum length of dimension
            """
            dim = self.dimensions[name] = PseudoNetCDFDimension(self, name, length)
            return dim

        def createVariable(self, name, type, dimensions, **properties):
            """
            name - string
            type - numpy dtype code (e.g., 'f', 'i', 'd')
            dimensions - tuple of dimension keys that can be
                         found in objects' dimensions dictionary
            """
            var = self.variables[name] = PseudoNetCDFVariable(self,name,type,dimensions, **properties)
            return var

        def close(self):
            """
            Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
            """
            pass

        def ncattrs(self):
            return self._ncattrs

        sync=close
        flush=close

    class PseudoNetCDFVariable(ndarray):
        """
        PseudoNetCDFVariable presents the Scientific.IO.NetCDF.NetCDFVariable interface,
        but unlike that type, provides a contructor for variables that could be used
        without adding it to the parent file
        """
        def __setattr__(self, k, v):
            """
            Set attributes (aka properties) and identify user-defined attributes.
            """
            if not hasattr(self, k) and k[:1] != '_':
                self._ncattrs += (k,)
            ndarray.__setattr__(self, k, v)
        def ncattrs(self):
            """
            Returns a tuple of attributes that have been user defined
            """
        
            return self._ncattrs
        def __new__(subtype,parent,name,typecode,dimensions,**kwds):
            """
            Creates a variable using the dimensions as defined in
            the parent object

            parent: an object with a dimensions variable
            name: name for variable
            typecode: numpy style typecode
            dimensions: a typle of dimension names to be used from
                        parrent
            kwds: Dictionary of keywords to be added as properties
                  to the variable.  **The keyword 'values' is a special
                  case that will be used as the starting values of
                  the array

            """
            if 'values' in kwds.keys():
                result=kwds.pop('values')
            else:
                shape=[]
                for d in dimensions:
                    dim = parent.dimensions[d]

                    # Adding support for netCDF3 dimension objects
                    if not isinstance(dim, int):
                        dim = len(dim)
                    shape.append(dim)

                result=zeros(shape,typecode)

            result=result[...].view(subtype)

            if hasattr(result, '__dict__'):
                result.__dict__['typecode'] = lambda: typecode
                result.__dict__['dimensions'] = tuple(dimensions)
            else:
                result.__dict__ = {
                    'typecode': lambda: typecode,
                    'dimensions': tuple(dimensions)
                }

    #        object.__setattr__(result, '_ncattrs', ())

            for k,v in kwds.iteritems():
                setattr(result,k,v)
            return result

        def __array_finalize__(self, obj):
            """
            finalization involves propagating features through opertations.
            """
            assert(hasattr(self, '_ncattrs') == False)
            self._ncattrs = ()
            if obj is None: return
            if hasattr(obj, '_ncattrs'):
                for k in obj._ncattrs:
                    setattr(self, k, getattr(obj, k))
            if not hasattr(self, 'dimensions'):
                if hasattr(obj, 'dimensions'):
                    setattr(self, 'dimensions', obj.dimensions)
                    self._ncattrs = self._ncattrs[:-1]
    
        def getValue(self):
            """
            Return scalar value
            """
            return self.item()

        def assignValue(self,value):
            """
            assign value to scalar variable
            """
            self.itemset(value)


# These variables define the binary format of the header blocks
# and are only for internal
_general_header_type = dtype('>i4, S40, >i4, >i4, S80, >i4')
_datablock_header_type = dtype('>i4, S20, 2>f4, >i4, >i4, >i4, >i4, S40, >i4, S40, >f8, >f8, S40, 3>i4, 3>i4, >i4, >i4')
_first_header_size = _general_header_type.itemsize + _datablock_header_type.itemsize

class defaultdictfromkey(defaultdict):
    """
    defaultdictfromkey dynamically produces dictionary items
    for keys, using the default_factor function called
    with the key as an argument
    """

    def __missing__(self, key):
        """
        __missing__(key) # Called by __getitem__ for missing key; pseudo-code:
        if self.default_factory is None: raise KeyError((key,))
        self[key] = value = self.default_factory()
        return value
        """
        return self.default_factory(key)

class defaultdictfromthesekeys(defaultdict):
    """
    defaultdictfromthesekeys dynamically produces dictionary items
    for known keys, using the default_factor function called
    with the key as an argument
    """
    def __init__(self, keys, default_factory = None):
        """
        keys - iterable of keys that default_factory can produce
        default_factory - function that takes a key as an argument
                          to create a dictionary item
        """
        self._keys = set([k for k in keys])
        defaultdict.__init__(self, default_factory)
    
    def __iter__(self):
        """
        Just like dictionary, but iterates through pre-defined keys
        """
        for i in self._keys:
            yield i
            
    def iterkeys(self):
        """
        Just like dictionary, but iterates through pre-defined keys
        """
        for i in self._keys:
            yield i
    
    def itervalues(self):
        """
        Just like dictionary, but iterates through pre-defined key values
        """
        for k in self.iterkeys():
            yield self[k]
    
    def iteritems(self):
        """
        Just like dictionary, but iterates through pre-defined keys and their calculated values
        """
        for k in self.iterkeys():
            yield (k, self[k])
    
    def keys(self):
        """
        Just like dictionary, but iterates through pre-defined keys
        """
        return [k for k in self]
    
    def __setitem__(self, key, value):
        """
        Add value with key and add to pre-defined keys
        """
        val = defaultdict.__setitem__(self, key, value)
        self._keys.add(key)
        return val
    
    def __delitem__(self, key):
        """
        Delete value with key and remove pre-defined key
        """
        self._keys.discard(key)
        return defaultdict.__delitem__(self, key)
    
    def pop(self, key):
        """
        Pop value with key and remove pre-defined key
        """
        val = defaultdict.pop(self, key)
        self._keys.discard(key)
        return val
        
    def __missing__(self, key):
        """
        __missing__(key) # Called by __getitem__ for missing key; pseudo-code:
        if self.default_factory is None: raise KeyError((key,))
        self[key] = value = self.default_factory()
        return value
        """
        if key in self._keys:
            return self.default_factory(key)
        else:
            raise KeyError("%s not found" % (key, ))

class defaultpseudonetcdfvariable(defaultdictfromthesekeys):
    """
    Overwrites __repr__ function to show variables
    """
    def __repr__(self):
        out = "{"
        for k in self._keys:
            out += "\n  '%s': PseudoNetCDFVariable(...)" % k
        out += "\n}"
        return out

    def __str__(self):
        return self.__repr__()

class _diag_group(PseudoNetCDFFile):
    """
    This object acts as a PseudoNetCDF file that gets data from the parent object.
    """
    def __init__(self, parent, groupname, groupvariables):
        """
        parent - a PseudoNetCDFFile
        groupname - a string describing the group
        """
        template = '%s_%%s' % groupname
        def getvar(key):
            try:
                return parent.variables[template % key]
            except (KeyError, ValueError):
                return parent.variables[key]
        self._parent = parent
        if 'BXHGHT-$_BXHEIGHT' not in groupvariables:
            mymetakeys = [k for k in metakeys if k != 'VOL']
        else:
            mymetakeys = metakeys
        self.variables = defaultpseudonetcdfvariable(list(groupvariables) + mymetakeys, getvar)
        
    
    def __getattr__(self, key):
        try:
            return object.__getattr__(self, key)
        except AttributeError, e:
            if key != 'groups':
                return getattr(self._parent, key)
            else:
                raise e
    
# This class is designed to operate like a dictionary, but
# dynamically create variables to return to the user
class _tracer_lookup(defaultpseudonetcdfvariable):
    """
    _tracer_lookup: finds geos_chem tracer indices from names and returns
                    netcdf like variable
    """
    def __init__(self, parent, datamap, tracerinfo, diaginfo, keys, noscale = False):
        """
        parent: NetCDF-like object to serve dimensions
        datamap: array of pre-dimensioned and datatyped values dim(tstep, i, j, k)
        tracerinfo: dictionary of tracer data keyed by ordinal
        keys: list of keys to serve
        """
        self._noscale = noscale
        self._tracer_data = tracerinfo
        self._diag_data = diaginfo
        self._memmap = datamap
        self._parent = parent
        self._special_keys = set(metakeys)
        self._keys = set(keys).union(self._special_keys)
        if 'BXHGHT-$_BXHEIGHT' not in keys:
            self._keys.discard('VOL')
        self._example_key = keys[0]
        
    def __missing__(self, key):
        if key in ('latitude', 'latitude_bounds'):
            yres = self._parent.modelres[1]
            if self._parent.halfpolar == 1:
                data = concatenate([[-90.], arange(-90. + yres / 2., 90., yres), [90.]])
            else:
                data = arange(-90, 90 + yres, yres)
            
            dims = ('latitude',)
            dtype = 'f'
            kwds = dict(standard_name = "latitude", long_name = "latitude", units = "degrees_north", base_units = "degrees_north", axis = "Y")
            if key == 'latitude':
                data = data[:-1] + diff(data) / 2.
                kwds['bounds'] = 'latitude_bounds'
            else:
                dims += ('nv',)
                data = data.repeat(2,0)[1:-1].reshape(-1, 2)
            example = self[self._example_key]
            sj = getattr(example, 'STARTJ', 0)
            data = data[sj:sj + example.shape[2]]
        elif key in ('longitude', 'longitude_bounds'):
            xres = self._parent.modelres[0]
            i = arange(0, 360 + xres, xres)
            data = i - (180 + xres / 2. * self._parent.center180)
            dims = ('longitude',)
            dtype = 'f'
            kwds = dict(standard_name = "longitude", long_name = "longitude", units = "degrees_east", base_units = "degrees_east", axis = "X")
            if key == 'longitude':
                data = data[:-1] + diff(data) / 2.
                kwds['bounds'] = 'longitude_bounds'
            else:
                dims += ('nv',)
                data = data.repeat(2,0)[1:-1].reshape(-1, 2)
            example = self[self._example_key]
            si = getattr(example, 'STARTI', 0)
            data = data[si:si + example.shape[3]]
        elif key == 'AREA':
           lon = self['longitude']
           xres = self._parent.modelres[0]
           nlon = 360. / xres
           latb = self['latitude_bounds']
           Re = self['crs'].semi_major_axis
           latb = append(latb[:, 0], latb[-1, 1])
           latr = pi / 180. * latb
           data = 2. * pi * Re * Re / (nlon) * ( sin( latr[1:] ) - sin( latr[:-1] ) )
           data = data[:, None].repeat(lon.size, 1)
           kwds = dict(units = 'm**2', base_units = 'm**2', grid_mapping = "crs")
           dtype = 'f'
           dims = ('latitude', 'longitude')
        elif key == 'VOL':
           try:
               bxhght = self['BXHGHT-$_BXHEIGHT']
               area = self['AREA']
           except KeyError:
               raise KeyError('Volume is only available if BXHGHT-$_BXHEIGHT was output')
           data = area[None,None] * bxhght
           kwds = dict(units = 'm**3', base_units = 'm**3', grid_mapping = "crs")
           dtype = 'f'
           dims = ('time', 'layer', 'latitude', 'longitude')
           if len(['layer' in dk_ for dk_ in self._parent.dimensions]) > 1:
               dims = ('time', 'layer%d' % data.shape[1], 'latitude', 'longitude')
        elif key == 'crs':
          dims = ()
          kwds = dict(grid_mapping_name = "latitude_longitude",
                      semi_major_axis = 6375000.0,
                      inverse_flattening = 0)
          dtype = 'i'
          data = zeros(1, dtype = dtype)
        elif key in ('time', 'time_bounds'):
            tmp_key = self._example_key
            data = np.array([self['tau0'], self['tau1']]).T
            dims = ('time', 'nv')
            if key == 'time':
                data = data.mean(1)
                dims = ('time',)
            
            dtype = 'd'
            kwds = dict(units = 'hours since 1985-01-01 00:00:00 UTC', base_units = 'hours since 1985-01-01 00:00:00 UTC', standard_name = key, long_name = key, var_desc = key)
            if key == 'time':
                kwds['bounds'] = 'time_bounds'
        elif key[:5] == 'layer':
            if self._parent.vertgrid in ('GEOS-5-REDUCED', 'MERRA-REDUCED'):
                if key[5:] in ('48', '_edges'):
                    data = np.array([1013.25, 998.051, 982.765, 967.48, 952.195, 936.911, 921.626, 906.342, 891.059, 875.776, 860.493, 845.211, 829.929, 809.556, 784.088, 758.621, 733.16, 707.699, 682.239, 644.054, 605.88, 567.706, 529.55, 491.401, 453.269, 415.155, 377.07, 339.005, 288.927, 245.246, 208.244, 176.93, 150.393, 127.837, 108.663, 92.366, 78.512, 56.388, 40.175, 28.368, 19.792, 9.293, 4.077, 1.651, 0.617, 0.211, 0.066, 0.01], dtype = 'f')
                else:
                    data = np.array([1005.65, 990.408, 975.122, 959.837, 944.553, 929.268, 913.984, 898.701, 883.418, 868.135, 852.852, 837.57, 819.743, 796.822, 771.354, 745.89, 720.429, 694.969, 663.146, 624.967, 586.793, 548.628, 510.475, 472.335, 434.212, 396.112, 358.038, 313.966, 267.087, 226.745, 192.587, 163.661, 139.115, 118.25, 100.514, 85.439, 67.45, 48.282, 34.272, 24.08, 14.542, 6.685, 2.864, 1.134, 0.414, 0.139, 0.038], dtype = 'f')
                    
            elif self._parent.vertgrid in ('GEOS-5-NATIVE', 'MERRA-NATIVE'):
                if key[5:] in ('73', '_edges'):
                    data = np.array([1013.25, 998.051, 982.765, 967.48, 952.195, 936.911, 921.626, 906.342, 891.059, 875.776, 860.493, 845.211, 829.929, 809.556, 784.088, 758.621, 733.16, 707.699, 682.239, 644.054, 605.88, 567.706, 529.55, 491.401, 453.269, 415.155, 377.07, 339.005, 288.927, 245.246, 208.244, 176.93, 150.393, 127.837, 108.663, 92.366, 78.512, 66.603, 56.388, 47.644, 40.175, 33.81, 28.368, 23.73, 19.792, 16.457, 13.643, 11.277, 9.293, 7.62, 6.217, 5.047, 4.077, 3.276, 2.62, 2.085, 1.651, 1.301, 1.019, 0.795, 0.617, 0.476, 0.365, 0.279, 0.211, 0.159, 0.12, 0.089, 0.066, 0.048, 0.033, 0.02, 0.01], dtype = 'f')
                else:
                    data = np.array([1005.65, 990.408, 975.122, 959.837, 944.553, 929.268, 913.984, 898.701, 883.418, 868.135, 852.852, 837.57, 819.743, 796.822, 771.354, 745.89, 720.429, 694.969, 663.146, 624.967, 586.793, 548.628, 510.475, 472.335, 434.212, 396.112, 358.038, 313.966, 267.087, 226.745, 192.587, 163.661, 139.115, 118.25, 100.514, 85.439, 72.558, 61.496, 52.016, 43.91, 36.993, 31.089, 26.049, 21.761, 18.124, 15.05, 12.46, 10.285, 8.456, 6.918, 5.632, 4.562, 3.677, 2.948, 2.353, 1.868, 1.476, 1.16, 0.907, 0.706, 0.546, 0.42, 0.322, 0.245, 0.185, 0.14, 0.105, 0.078, 0.057, 0.04, 0.026, 0.015], dtype = 'f')
            elif self._parent.vertgrid == 'GEOS-4-REDUCED':
                data = np.array([1005.706, 983.328, 941.644, 877.974, 797.026, 704.433, 606.413, 514.754, 436.898, 370.793, 314.683, 267.087, 226.745, 192.587, 163.661, 139.115, 118.25, 100.515, 85.439, 67.45, 48.282, 34.272, 24.08, 14.542, 6.685, 2.864, 1.134, 0.414, 0.139, 0.038], dtype = 'f')
            elif self._parent.vertgrid == 'GEOS-4-NATIVE':
                data = np.array([1005.706, 983.328, 941.644, 877.974, 797.026, 704.433, 606.413, 514.754, 436.898, 370.793, 314.683, 267.087, 226.745, 192.587, 163.661, 139.115, 118.25, 100.515, 85.439, 72.558, 61.496, 52.016, 43.91, 36.993, 31.089, 26.049, 21.761, 18.124, 15.05, 12.46, 10.285, 8.456, 6.918, 5.632, 4.562, 3.677, 2.948, 2.353, 1.868, 1.476, 1.16, 0.907, 0.706, 0.546, 0.42, 0.322, 0.245, 0.185, 0.14, 0.105, 0.078, 0.057, 0.04, 0.026, 0.015], dtype = 'f')
            else:
                warn('Assuming GEOS-5 reduced vertical coordinate')
                data = np.array([1005.65, 990.408, 975.122, 959.837, 944.553, 929.268, 913.984, 898.701, 883.418, 868.135, 852.852, 837.57, 819.743, 796.822, 771.354, 745.89, 720.429, 694.969, 663.146, 624.967, 586.793, 548.628, 510.475, 472.335, 434.212, 396.112, 358.038, 313.966, 267.087, 226.745, 192.587, 163.661, 139.115, 118.25, 100.514, 85.439,67.45, 48.282, 34.272, 24.08, 14.542, 6.685, 2.864, 1.134, 0.414, 0.139, 0.038], dtype = 'f')
            
            data = data[:len(self._parent.dimensions[key])]
            dims = (key,)
            dtype = 'f'
            kwds = dict(units = 'model layer', base_units = 'model layer', standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate', long_name = key, var_desc = key, axis = "Z")
        elif key == 'tau0':
            tmp_key = self._example_key
            data = self._memmap[tmp_key]['header']['f10']
            dims = ('time',)
            dtype = 'd'
            kwds = dict(units = 'hours since 1985-01-01 00:00:00 UTC', base_units = 'hours since 1985-01-01 00:00:00 UTC', standard_name = key, long_name = key, var_desc = key)
        elif key == 'tau1':
            tmp_key = self._example_key
            data = self._memmap[tmp_key]['header']['f11']
            dims = ('time',)
            dtype = 'd'
            kwds = dict(units = 'hours since 1985-01-01 00:00:00 UTC', base_units = 'hours since 1985-01-01 00:00:00 UTC', standard_name = key, long_name = key, var_desc = key)
        else:
            dtype = 'f'
            header = self._memmap[key]['header'][0]
            sl, sj, si = header['f14'][::-1] - 1
            group = header['f7'].strip()
            offset = self._diag_data.get(group, {}).get('offset', 0)
            ord = header['f8'] + offset
            base_units = header['f9']
            scale = self._tracer_data[ord]['SCALE']
            molwt = self._tracer_data[ord]['MOLWT']
            carbon = self._tracer_data[ord]['C']
            units = self._tracer_data[ord]['UNIT']
            tmp_data = self._memmap[key]['data']
            dims = ('time', 'layer', 'latitude', 'longitude')
            if len(['layer' in dk_ for dk_ in self._parent.dimensions]) > 1:
                dims = ('time', 'layer%d' % tmp_data.dtype['f1'].shape[0], 'latitude', 'longitude')
            kwds = dict(scale = scale, kgpermole = molwt, carbon = carbon, units = units, base_units = base_units, standard_name = key, long_name = key, var_desc = key, coordinates = ' '.join(dims), grid_mapping = "crs")
                
            assert((tmp_data['f0'] == tmp_data['f2']).all())
            if self._noscale:
                if scale != 1.:
                    warn("Not scaling variables; good for writing")
                data = tmp_data['f1']
            else:
                data = tmp_data['f1'] * scale
            if any([sl != 0, sj != 0, si != 0]):
                nl, nj, ni = header['f13'][::-1]
                #import pdb; pdb.set_trace()
                #tmp_data = zeros((data.shape[0], self._parent.dimensions['layer'], self._parent.dimensions['latitude'], self._parent.dimensions['longitude']), dtype = data.dtype)
                #el, ej, ei = data.shape[1:]
                #el += sl
                #ej += sj
                #ei += si
                #tmp_data[:, sl:el, sj:ej, si:ei] = data[:]
                #data = tmp_data
                kwds['STARTI'] = si 
                kwds['STARTJ'] = sj
                kwds['STARTK'] = sl
        return PseudoNetCDFVariable(self._parent, key, dtype, dims, values = data, **kwds)

coordkeys = 'time latitude longitude layer latitude_bounds longitude_bounds crs'.split()            
metakeys = ['VOL', 'AREA', 'tau0', 'tau1', 'time_bounds'] + coordkeys


class bpch(PseudoNetCDFFile):
    """
    NetCDF-like class to interface with GEOS-Chem binary punch files
    
    f = bpch(path_to_binary_file)
    dim = f.dimensions[dkey] # e.g., dkey = 'longitude'
    
    # There are two ways to get variables.  Directly from
    # the file using the long name
    var = f.variables[vkey] # e.g., vkey = 'IJ-AVG-$_NOx'
    
    # Or through a group using the short name
    g = f.groups[gkey] # e.g., gkey = 'IJ-AVG-$'
    var = g.variables[vkey] # e.g., vkey = 'NOx'

    # The variable returned is the same either way    
    print f.dimensions
    print var.unit
    print var.dimensions
    print var.shape
    
    """

    def __init__(self, bpch_path, tracerinfo = None, diaginfo = None, mode = 'r', timeslice = slice(None), noscale = False, vertgrid = 'GEOS-5-REDUCED'):
        """
        bpch_path: path to binary punch file
        tracerinfo: path to ascii file with tracer definitions
        diaginfo: path to ascii file with diagnostic group definitions
        mode : {'r+', 'r', 'w+', 'c'}, optional
         |      The file is opened in this mode:
         |  
         |      +------+-------------------------------------------------------------+
         |      | 'r'  | Open existing file for reading only.                        |
         |      +------+-------------------------------------------------------------+
         |      | 'r+' | Open existing file for reading and writing.                 |
         |      +------+-------------------------------------------------------------+
         |      | 'w+' | Create or overwrite existing file for reading and writing.  |
         |      +------+-------------------------------------------------------------+
         |      | 'c'  | Copy-on-write: assignments affect data in memory, but       |
         |      |      | changes are not saved to disk.  The file on disk is         |
         |      |      | read-only.                                                  |
         |      +------+-------------------------------------------------------------+        
         timeslice: If the file is larger than 2GB, timeslice provides a way to subset results.
                    The subset requested depends on the data type of timeslice:
                        - int: return the a part of the file if it was broken into 2GB chunks (0..N-1)
                        - slice: return the times that correspond to that slice (i.e., range(ntimes)[timeslice])
                        - list/tuple/set: return specified times where each time is in the set (0..N-1)
         noscale: Do not apply scaling factors
         vertgrid: vertical coordinate system (options: 'GEOS-5-REDUCED', 'GEOS-5-NATIVE', 'MERRA-REDUCED', 'MERRA-NATIVE', 'GEOS-4-REDUCED', 'GEOS-4-NATIVE' -- default 'GEOS-5-REDUCED')
        """
        self._ncattrs = () 
        self._noscale = noscale
        # Read binary data for general header and first datablock header
        header_block = fromfile(bpch_path, dtype = 'bool', count = _first_header_size)
        
        # combine data for convenience
        header = tuple(header_block[:_general_header_type.itemsize].view(_general_header_type)[0]) + \
                 tuple(header_block[_general_header_type.itemsize:].view(_datablock_header_type)[0])
        
        # Verify that all Fortran unformatted buffers match 
        try:
            assert(header[0] == header[2])
            assert(header[3] == header[5])
        except AssertionError:
            raise ValueError("BPCH Files fails header check")
        
        # Assign data from header to global attributes
        self.ftype = header[1]
        self.toptitle = header[4]
        self.modelname, self.modelres, self.halfpolar, self.center180 = header[7:11]
        dummy, dummy, dummy, self.start_tau0, self.start_tau1, dummy, dim, dummy, dummy = header[13:-1]
        for dk, dv in zip('longitude latitude layer'.split(), dim):
            self.createDimension(dk, dv)
        self.createDimension('nv', 2)
        tracerinfo = tracerinfo or os.path.join(os.path.dirname(bpch_path), 'tracerinfo.dat')
        if isinstance(tracerinfo, (str, unicode)):
            if not os.path.exists(tracerinfo) and tracerinfo != ' ':
                tracerinfo = 'tracerinfo.dat'
            if os.path.exists(tracerinfo):
                if os.path.isdir(tracerinfo): tracerinfo = os.path.join(tracerinfo, 'tracerinfo.dat')
        
            tracerinfo = file(tracerinfo)
        if isinstance(tracerinfo, (file, StringIO)):
            tracer_data = dict([(int(l[52:61].strip()), dict(NAME = l[:8].strip(), FULLNAME = l[9:39].strip(), MOLWT = float(l[39:49]), C = int(l[49:52]), TRACER = int(l[52:61]), SCALE = float(l[61:71]), UNIT = l[72:].strip())) for l in tracerinfo.readlines() if l[0] not in ('#', ' ')])
            tracer_names = dict([(k, v['NAME']) for k, v in tracer_data.iteritems()])
        else:
            warn('Reading file without tracerinfo.dat means that names and scaling are unknown')
            tracer_data = defaultdict(lambda: dict(SCALE = 1., C = 1.))
            tracer_names = defaultdictfromkey(lambda key: key)
        
        diaginfo = diaginfo or os.path.join(os.path.dirname(bpch_path), 'diaginfo.dat')
        if isinstance(diaginfo, (str, unicode)):
            if not os.path.exists(diaginfo):
                diaginfo = 'diaginfo.dat'
            if os.path.exists(diaginfo):
                if os.path.isdir(diaginfo): diaginfo = os.path.join(diaginfo, 'diaginfo.dat')
            diaginfo = file(diaginfo)

        if isinstance(diaginfo, (file, StringIO)):
            diag_data = dict([(l[9:49].strip(), dict(offset = int(l[:8]), desc = l[50:].strip())) for l in diaginfo.read().strip().split('\n') if l[0] != '#'])
        else:
            warn('Reading file without diaginfo.dat loses descriptive information')
            diag_data = defaultdictfromkey(lambda key: dict(offset = 0, desc = key))
            
        if len(tracer_names) == 0 and not isinstance(tracer_names, defaultdictfromkey):
            raise IOError("Error parsing %s for Tracer data" % tracerinfo)
        file_size = os.stat(bpch_path).st_size
        offset = _general_header_type.itemsize
        data_types = []
        first_header = None
        keys = []
        self._groups = defaultdict(set)
        while first_header is None or \
              offset < file_size:
            header = memmap(bpch_path, offset = offset, shape = (1,), dtype = _datablock_header_type, mode = mode)[0]
            
            group = header[7].strip()
            tracer_number = header[8]
            unit = header[9].strip()
            if not isinstance(diag_data, defaultdictfromkey):
                goffset = diag_data.get(group, {}).get('offset', 0)
                try:
                    tracername = tracer_names[tracer_number + goffset]
                except:
                    # There are some cases like with the adjoint where the tracer does not have
                    # a tracerinfo.dat line.  In this case, the name matches the tracer with that
                    # number (no offset).  The scaling, however, is not intended to be used.
                    # The unit for adjoint, for instance, is unitless.
                    if tracer_number not in tracer_names:
                        tracername = str(tracer_number)
                        tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 0, UNIT = unit)
                    else:
                        tracername = tracer_names[tracer_number]
                        tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = tracer_data[tracer_number]['C'], UNIT = unit)

            else:
                warn('%s is not in diaginfo.dat; names and scaling cannot be resolved' % group)
                goffset = 0
                tracername = str(tracer_number)
                diag_data[group] = dict(offset = 0, desc = group)
                tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 1., UNIT = unit)

            self._groups[group].add(tracername)
            offset += _datablock_header_type.itemsize + header[-2]

            if first_header is None:
                first_header = header
            elif (header[7], header[8]) == (first_header[7], first_header[8]) or offset == file_size:
                if offset == file_size:
                    dim = header[13][::-1]
                    start = header[14][::-1]
                    data_type = dtype('>i4, %s>f4, >i4' % str(tuple(dim[:])))
                    assert(data_type.itemsize == header[-2])
                    data_types.append(data_type)
                    keys.append('%s_%s' % (group, tracername))

                # Repeating tracer indicates end of timestep
                time_type = dtype([(k, dtype([('header', _datablock_header_type), ('data', d)])) for k, d in zip(keys, data_types)])
                field_shapes = set([v[0].fields['data'][0].fields['f1'][0].shape for k, v in time_type.fields.iteritems()])
                field_levs = set([s_[0] for s_ in field_shapes])
                field_rows = set([s_[1] for s_ in field_shapes])
                field_cols = set([s_[2] for s_ in field_shapes])
                field_levs = list(field_levs)
                field_levs.sort()
                for fl in field_levs:
                    self.createDimension('layer%d' % fl, fl)

                itemcount = ((float(os.path.getsize(bpch_path)) - _general_header_type.itemsize) / time_type.itemsize)
                if (itemcount % 1) != 0:
                    warn("Cannot read whole file; assuming partial time block is at the end; skipping partial time record")
                    itemcount = np.floor(itemcount)
                break
            dim = header[13][::-1]
            start = header[14][::-1]
            data_type = dtype('>i4, %s>f4, >i4' % str(tuple(dim[:])))
            assert(data_type.itemsize == header[-2])
            data_types.append(data_type)
            keys.append('%s_%s' % (group, tracername))

        # load all data blocks  
        try:
            datamap = memmap(bpch_path, dtype = time_type, offset = _general_header_type.itemsize, mode = mode, shape = (itemcount,))
            if not timeslice is None:
                datamap = datamap[timeslice]
        except OverflowError:
            hdrsize = _general_header_type.itemsize
            items = (2*1024**3-hdrsize) // time_type.itemsize
            if timeslice != slice(None):
                filesize = os.stat(bpch_path).st_size
                datasize = (filesize - hdrsize)
                all_times = range(datasize / time_type.itemsize)
                if isinstance(timeslice, int):
                    timeslice = slice(items * timeslice, items * (timeslice + 1))
                if isinstance(timeslice, (list, tuple, set)):
                    times = timeslice
                else:
                    times = all_times[timeslice]

                outpath = bpch_path + '.tmp.part'
                mint = times[0]
                maxt = times[-1]
                nt = maxt - mint + 1
                
                if nt > items:
                    warn('Requested %d items; only returning %d items due to 2GB limitation' % (nt, items))
                    times = times[:items]

                outfile = file(outpath, 'w')
                infile = file(bpch_path, 'r')
                hdr = infile.read(hdrsize)
                outfile.write(hdr)
                for t in all_times:
                    if t in times:
                        outfile.write(infile.read(time_type.itemsize))
                    else:
                        infile.seek(time_type.itemsize, 1)
                outfile.close()
                #print mint, maxt, nt, nt * time_type.itemsize
                #cmd = 'dd if=%s ibs=%d skip=1 obs=%d | dd of=%s bs=%d skip=%d count=%d' % (bpch_path, hdrsize, time_type.itemsize, outpath, time_type.itemsize, mint, nt)
                #print cmd
                #os.system(cmd)
                datamap = memmap(outpath, dtype = time_type, mode = mode, offset = hdrsize)
            else:
                datamap = memmap(bpch_path, dtype = time_type, shape = (items,), offset = _general_header_type.itemsize, mode = mode)
                warn('Returning only the first 2GB of data')
        for k in datamap.dtype.names:
          gn = datamap[k]['header']['f7']
          tid = datamap[k]['header']['f8']
          assert((tid[0] == tid).all())
          assert((gn[0] == gn).all())
        # Create variables and dimensions
        self.vertgrid = vertgrid
        layerns = set([datamap[0][k]['header']['f13'][-1] for k in datamap.dtype.names])
        layerkeys = ['layer_edges'] + ['layer%d' % l for l in layerns]
        keys.extend(layerkeys)
        if self.vertgrid in ('GEOS-5-REDUCED', 'MERRA-REDUCED'):
            self.createDimension('layer', 47)
            self.createDimension('layer_edges', 48)
            # Ap [hPa] for 47 levels (48 edges)
            self.Ap = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
                             1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
                             4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
                             7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
                             1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
                             1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
                             2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
                             2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
                             1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
                             7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, 
                             1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00, 
                             6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])

            # Bp [unitless] for 47 levels (48 edges)
            self.Bp = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
                             9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
                             8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
                             7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
                             6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
                             4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
                             2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
                             6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

        elif self.vertgrid in ('GEOS-5-NATIVE', 'MERRA-NATIVE'):
            self.createDimension('layer', 72)
            self.createDimension('layer_edges', 73)
            # Ap [hPa] for 72 levels (73 edges)
            self.Ap = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
                     1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
                     4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
                     7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
                     1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
                     1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
                     2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
                     2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
                     1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
                     7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01,
                     4.017541e+01, 3.381001e+01, 2.836781e+01, 2.373041e+01,
                     1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
                     9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00,
                     4.076571e+00, 3.276431e+00, 2.620211e+00, 2.084970e+00,
                     1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
                     6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01,
                     2.113490e-01, 1.594950e-01, 1.197030e-01, 8.934502e-02,
                     6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
                     1.000000e-02])

            # Bp [unitless] for 72 levels (73 edges)
            self.Bp = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
                     9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
                     8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
                     7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
                     6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
                     4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
                     2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
                     6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00])
        elif self.vertgrid == 'GEOS-4-REDUCED':
            self.createDimension('layer', 30)
            self.createDimension('layer_edges', 31)
            # Ap [hPa] for 30 levels (31 edges)
            self.AP = np.array([0.000000e0,   0.000000e0,  12.704939e0,  35.465965e0, 
                        66.098427e0, 101.671654e0, 138.744400e0, 173.403183e0, 
                       198.737839e0, 215.417526e0, 223.884689e0, 224.362869e0, 
                       216.864929e0, 201.192093e0, 176.929993e0, 150.393005e0, 
                       127.837006e0, 108.663429e0,  92.365662e0,  78.512299e0, 
                        56.387939e0,  40.175419e0,  28.367815e0,  19.791553e0, 
                         9.292943e0,   4.076567e0,   1.650792e0,   0.616779e0, 
                         0.211349e0,   0.066000e0,   0.010000e0])

            # Bp [unitless] for 30 levels (31 edges)
            self.Bp = np.array([1.000000e0,   0.985110e0,   0.943290e0,   0.867830e0, 
                         0.764920e0,   0.642710e0,   0.510460e0,   0.378440e0, 
                         0.270330e0,   0.183300e0,   0.115030e0,   0.063720e0, 
                         0.028010e0,   0.006960e0,   0.000000e0,   0.000000e0, 
                         0.000000e0,   0.000000e0,   0.000000e0,   0.000000e0, 
                         0.000000e0,   0.000000e0,   0.000000e0,   0.000000e0, 
                         0.000000e0,   0.000000e0,   0.000000e0,   0.000000e0, 
                         0.000000e0,   0.000000e0,   0.000000e0])
        elif self.vertgrid == 'GEOS-4-NATIVE':
            self.createDimension('layer', 55)
            self.createDimension('layer_edges', 56)
            # AP [hPa] for 55 levels (56 edges)
            self.Ap = np.array([0.000000e0,   0.000000e0,  12.704939e0,  35.465965e0, 
                       66.098427e0, 101.671654e0, 138.744400e0, 173.403183e0,
                      198.737839e0, 215.417526e0, 223.884689e0, 224.362869e0,
                      216.864929e0, 201.192093e0, 176.929993e0, 150.393005e0,
                      127.837006e0, 108.663429e0,  92.365662e0,  78.512299e0, 
                       66.603378e0,  56.387939e0,  47.643932e0,  40.175419e0, 
                       33.809956e0,  28.367815e0,  23.730362e0,  19.791553e0, 
                       16.457071e0,  13.643393e0,  11.276889e0,   9.292943e0,
                        7.619839e0,   6.216800e0,   5.046805e0,   4.076567e0, 
                        3.276433e0,   2.620212e0,   2.084972e0,   1.650792e0,
                        1.300508e0,   1.019442e0,   0.795134e0,   0.616779e0, 
                        0.475806e0,   0.365041e0,   0.278526e0,   0.211349e0, 
                        0.159495e0,   0.119703e0,   0.089345e0,   0.066000e0, 
                        0.047585e0,   0.032700e0,   0.020000e0,   0.010000e0])

            # BP [unitless] for 55 levels (56 edges)
            self.Bp = np.array([1.000000e0,  0.985110e0,   0.943290e0,   0.867830e0,
                         0.764920e0,  0.642710e0,   0.510460e0,   0.378440e0,
                         0.270330e0,  0.183300e0,   0.115030e0,   0.063720e0,
                         0.028010e0,  0.006960e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0,
                         0.000000e0,  0.000000e0,   0.000000e0,   0.000000e0])
        else:
            warn('Assuming GEOS-5 reduced vertical coordinate')
            self.createDimension('layer', 47)
            self.createDimension('layer_edges', 48)
            self.Ap = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
                             1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
                             4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
                             7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
                             1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
                             1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
                             2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
                             2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
                             1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
                             7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, 
                             1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00, 
                             6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])

            # Bp [unitless] for 47 levels (48 edges)
            self.Bp = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
                             9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
                             8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
                             7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
                             6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
                             4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
                             2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
                             6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

        self.variables = _tracer_lookup(parent = self, datamap = datamap, tracerinfo = tracer_data, diaginfo = diag_data, keys = keys, noscale = self._noscale)
        del datamap
        tdim = self.createDimension('time', self.variables['tau0'].shape[0])
        tdim.setunlimited(True)
        self.groups = dict([(k, _diag_group(self, k, v)) for k, v in self._groups.iteritems()])
        for grp in self.groups.values():
            for dk, d in self.dimensions.iteritems():
                dmn = grp.createDimension(dk, len(d))
                dmn.setunlimited(d.isunlimited())
        self.Conventions = "CF-1.6"


    def __repr__(self):
        return PseudoNetCDFFile.__repr__(self) + str(self.variables)

def tileplot(f, toplot, vmin = None, vmax = None, xmin = None, xmax = None, ymin = None, ymax = None, title = '', unit = '', maptype = 0, log = False, cmap = None):
    # Example: spatial plotting
    from pylab import figure, colorbar, axis
    from matplotlib.colors import Normalize, LogNorm
    has_map = False
    lat = np.append(f.variables['latitude_bounds'][0, 0], f.variables['latitude_bounds'][:, 1])
    lon = np.append(f.variables['longitude_bounds'][0, 0], f.variables['longitude_bounds'][:, 1])
    parallels = arange(lat.min(),lat.max() + 15,15)
    meridians = arange(lon.min(), lon.max() + 30,30)
    fig = figure(figsize = (9,4))
    ax = fig.add_axes([.1, .1, .9, .8])
    m = ax
    if maptype == 0:
        try:
            from mpl_toolkits.basemap import Basemap
    
            # I'm actually not sure what the right projection for this 
            # data is.  So the next few lines might need to change.
            m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-180 - f.modelres[1] * f.center180 / 2.,\
                    urcrnrlon=180 - f.modelres[1] * f.center180 / 2.,\
                    resolution='c')
            has_map = True
            x, y = np.meshgrid(*m(lon, lat))
        except Exception, e:
            warn(str(e))
            x, y = np.meshgrid(lon, lat)
    
        if has_map:
            m.drawcountries()
            m.drawstates()
            m.drawcoastlines()
            m.drawparallels(parallels)
            m.drawmeridians(meridians)
        x = lon
        y = lat
        ax.set_xticks(meridians[::2])
        ax.set_yticks(parallels[::2])
        ax.set_xlim(meridians.min(), meridians.max())
        ax.set_ylim(parallels.min(), parallels.max())
    elif maptype == 1:
        x = lat
        ax.set_xticks(parallels)
        y = np.arange(toplot.shape[0])
        if 'BXHGHT-$_BXHEIGHT' in f.variables.keys():
            y = np.cumsum(np.append(0, f.variables['BXHGHT-$_BXHEIGHT'].mean(0).mean(1).mean(1))) / 1000.
            y = y[:toplot.shape[0]]
            ax.set_ylabel('km agl')
        else:
            ax.set_ylabel('layer')
        ax.set_xlabel('degrees north')
    elif maptype == 2:
        x = lon
        y = np.arange(toplot.shape[0])
        ax.set_xticks(meridians)
        if 'BXHGHT-$_BXHEIGHT' in f.variables.keys():
            y = np.cumsum(np.append(0, f.variables['BXHGHT-$_BXHEIGHT'].mean(0).mean(1).mean(1))) / 1000.
            y = y[:toplot.shape[0]]
            ax.set_ylabel('km')
        else:
            ax.set_ylabel('layer')
        ax.set_xlabel('degrees east')
    elif maptype == 3:
        from datetime import datetime, timedelta
        x = np.append(f.variables['tau0'][0], f.variables['tau1'][:]) / 24.
        start_date = datetime(1985, 1, 1) + timedelta(hours = f.variables['tau0'][0])
        x = x - x.min()
        y = np.arange(toplot.shape[0])
        if 'BXHGHT-$_BXHEIGHT' in f.variables.keys():
            y = np.cumsum(np.append(0, f.variables['BXHGHT-$_BXHEIGHT'].mean(0).mean(1).mean(1))) / 1000.
            y = y[:toplot.shape[0]]
            ax.set_ylabel('km')
        else:
            ax.set_ylabel('layer')
        ax.set_xticks(np.linspace(x[0], x[-1], min(len(x), len(ax.get_xticks()))))
        ax.set_xlabel('days since %s' % start_date.strftime('%Y%m%dT%H%M%S'))
        toplot = toplot.T
    elif maptype == 4:
        from datetime import datetime, timedelta
        x = lat
        ax.set_xlabel('degrees north')
        ax.set_xticks(parallels)
        y = np.append(f.variables['tau0'][0], f.variables['tau1'][:]) / 24.
        start_date = datetime(1985, 1, 1) + timedelta(hours = f.variables['tau0'][0])
        y = y - y.min()
        ax.set_yticks(np.linspace(y[0], y[-1], min(len(y), len(ax.get_yticks()))))
        ax.set_ylabel('days since %s' % start_date.strftime('%Y%m%dT%H%M%S'))
    elif maptype == 5:
        from datetime import datetime, timedelta
        x = lon
        ax.set_xlabel('degrees east')
        ax.set_xticks(meridians)
        y = np.append(f.variables['tau0'][0], f.variables['tau1'][:]) / 24.
        start_date = datetime(1985, 1, 1) + timedelta(hours = f.variables['tau0'][0])
        y = y - y.min()
        ax.set_yticks(np.linspace(y[0], y[-1], min(len(y), len(ax.get_yticks()))))
        ax.set_ylabel('days since %s' % start_date.strftime('%Y%m%dT%H%M%S'))
    
    poly = m.pcolor(x, y, toplot, cmap = cmap)
    ax.collections[-1].set_norm((LogNorm if log else Normalize)(vmin = vmin, vmax = vmax))
    try:
        cb = colorbar(poly, ax = ax)
        cb.ax.set_xlabel(toplot.units.strip())
        if vmax is None:
            vmax = toplot.max()
        if vmin is None:
            vmin = np.ma.masked_values(toplot, 0).min()
        if log:
            if (np.log10(vmax) - np.log10(vmin)) < 4:
                ticks = np.logspace((np.log10(vmin)), (np.log10(vmax)), 10.)
                cb.set_ticks(ticks)
                cb.set_ticklabels(['%.1f' % x for x in ticks])
    except:
        pass
    axis('tight')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title)
    return fig

def getvar(path_to_test_file = '', group_key = None, var_key = None):
    # Example: file open and variable aquisition
    f = None
    while f is None:
        try:
            f = bpch(path_to_test_file)
        except Exception as e:
            print path_to_test_file, str(e)
            path_to_test_file = raw_input('Enter path to a valid GEOS-Chem file\n:')
    if group_key is None:
        group_names = f.groups.keys()
        if len(group_names) == 1:
            group_key, = group_names
        else:
            group_names.sort()
            while group_key not in group_names:
                group_key = raw_input('Enter a group name: %s\n:' % ', '.join(group_names))
    g = f.groups[group_key]
    
    if var_key is None:
        var_names = [k for k in g.variables.keys() if k not in metakeys]
        if len(var_names) == 1:
            var_key, = var_names
        else:
            var_names.sort()
            var_names = map(str, var_names)
            while var_key not in var_names:
                var_key = raw_input('Enter a variable name: %s\n:' % ', '.join(map(str, var_names)))

    try:
       var = g.variables[var_key]
    except:
       var = g.variables[int(var_key) if var_key.isdigit() else var_key]
    return f, group_key, var_key, var

def pad(nplots, option, option_key, default):
    nopts = len(option)
    if nplots != nopts:
        if nopts == 1:
            option = option * nplots
        else:
            if nopts != 0:
                warn("Padding %s with default layer = %s" % (option_key, str(default)))
            
            option = option + [default] * (nplots - nopts)
    return option

def reduce_dimold(var, eval_str, axis):
    if eval_str == 'animate':
        return var
    if eval_str in dir(np):
        var = getattr(np, eval_str)(var[:], axis = axis)[(slice(None),) * axis + (None,)]
    elif eval_str.isdigit():
        eval_idx = int(eval_str)
        var = np.take(var, [eval_idx], axis = axis)
    else:
        eval_idx = eval(eval_str)
        if isinstance(eval_idx, tuple) and len(eval_idx) == 3:
            var = getattr(np, eval_idx[0])(var[(slice(None),) * axis + (slice(eval_idx[1], eval_idx[2]),)], axis = axis)[(slice(None),) * axis + (None,)]
        else:
            var = var[(slice(None),) * axis + (eval_idx,)]
    return var

def slice_dim(f, slicedef, fuzzydim = True):
    slicedef = slicedef.split(',')
    slicedef = [slicedef[0]] + map(eval, slicedef[1:])
    if len(slicedef) == 2:
        slicedef.append(slicedef[-1] + 1)
    slicedef = (slicedef + [None,])[:4]
    dimkey, dmin, dmax, dstride = slicedef    
    if fuzzydim:
        partial_check = [key for key in f.dimensions if dimkey == key[:len(dimkey)] and key[len(dimkey):].isdigit()]
        for dimk in partial_check:
            f = slice_dim(f, '%s,%s,%s,%s' % (dimk, dmin, dmax, dstride))
        
    
    for varkey in f.variables.keys():
        var = f.variables[varkey]
        if dimkey not in var.dimensions:
            continue
        else:
            thisdimkey = dimkey
        axis = list(var.dimensions).index(dimkey)
        vout = var[:].swapaxes(0, axis)[dmin:dmax:dstride].swapaxes(0, axis)
        
        newlen = vout.shape[axis]
        f.createDimension(dimkey, newlen)
        f.variables[varkey] = vout
    return f
    
def reduce_dim(f, reducedef, fuzzydim = True):
    commacount = reducedef.count(',')
    if commacount == 3:
        dimkey, func, numweightkey, denweightkey = reducedef.split(',')
        numweight = f.variables[numweightkey]
        denweight = f.variables[denweightkey]
    elif commacount == 2:
        dimkey, func, numweightkey = reducedef.split(',')
        numweight = f.variables[numweightkey]
        denweightkey = None
    elif commacount == 1:
        dimkey, func = reducedef.split(',')
        numweightkey = None
        denweightkey = None

    if fuzzydim:
        partial_check = [key for key in f.dimensions if dimkey == key[:len(dimkey)] and key[len(dimkey):].isdigit()]
        for dimk in partial_check:
            if commacount == 1:
                f = reduce_dim(f, '%s,%s' % (dimk, func),)
            elif commacount == 2:
                f = reduce_dim(f, '%s,%s,%s' % (dimk, func, numweightkey),)
            elif commacount == 3:
                f = reduce_dim(f, '%s,%s,%s,%s' % (dimk, func, numweightkey, denweightkey),)
    
    f.createDimension(dimkey, 1)
    for varkey in f.variables.keys():
        var = f.variables[varkey]
        if dimkey not in var.dimensions:
            continue
        
        axis = list(var.dimensions).index(dimkey)
        
        if not varkey in metakeys:
            if numweightkey is None:
                vout = getattr(np, func)(var, axis = axis)[(slice(None),) * axis + (None,)]
            elif denweightkey is None:
                vout = getattr(np, func)(var * np.array(numweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)], axis = axis)[(slice(None),) * axis + (None,)]
                vout.units = vout.units.strip() + ' * ' + numweight.units.strip()
                if hasattr(vout, 'base_units'):
                    vout.base_units = vout.base_units.strip() + ' * ' + numweight.base_units.strip()
            else:
                vout = getattr(np, func)(var * np.array(numweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)], axis = axis)[(slice(None),) * axis + (None,)] / getattr(np, func)(np.array(denweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)], axis = axis)[(slice(None),) * axis + (None,)]
        else:
            if '_bounds' not in varkey:
                vout = getattr(np, func)(var, axis = axis)[(slice(None),) * axis + (None,)]
            else:
                vout = getattr(np, func)(var, axis = axis)[(slice(None),) * axis + (None,)]
                vout[0] = var[:].min(), var[:].max()
        f.variables[varkey] = vout
    return f

def getdiffpnc(f1, f2, percent = False):
    outf = getvarpnc(f1, None)
    f2 = getvarpnc(f2, None)
    for varkey in outf.variables.keys():
        if varkey not in coordkeys:
            outf.variables[varkey] = f1.variables[varkey] - f2.variables[varkey]
            if percent:
                outf.variables[varkey] /= f2.variables[varkey] * 100.
    return outf
    
def getvarpnc(f, varkeys):
    if varkeys is None:
        varkeys = list(set(f.variables.keys()).difference(coordkeys))
    else:
        newvarkeys = list(set(varkeys).intersection(f.variables.keys()))
        newvarkeys.sort()
        oldvarkeys = list(varkeys)
        oldvarkeys.sort()
        if newvarkeys != oldvarkeys:
            warn('Skipping %s' % ', '.join(set(oldvarkeys).difference(newvarkeys)))
        varkeys = newvarkeys

    outf = PseudoNetCDFFile()
    outf.createDimension('nv', 2)
    for propkey in f.ncattrs():
        setattr(outf, propkey, getattr(f, propkey))
    thiscoordkeys = [k for k in coordkeys]
    for varkey in varkeys:
        outf.DATAKEY = varkey
        try:
            var = eval(varkey, None, f.variables)
        except:
            var = f.variables[varkey]
        for dimk, dimv in zip(var.dimensions, var.shape):
            if dimk not in outf.dimensions:
                newdimv = outf.createDimension(dimk, dimv)
                if f.dimensions[dimk].isunlimited():
                    newdimv.setunlimited(True)
                coordkeys.append(dimk)
        for coordk in coordkeys:
            if coordk in f.dimensions and coordk not in outf.dimensions:
                newdimv = outf.createDimension(coordk, len(f.dimensions[coordk]))
                if f.dimensions[coordk].isunlimited():
                    newdimv.setunlimited(True)
    
        propd = dict([(k, getattr(var, k)) for k in var.ncattrs()])
        outf.createVariable(varkey, var.dtype.char, var.dimensions, values = var[:], **propd)
    for coordkey in coordkeys:
        coordvar = f.variables[coordkey]
        propd = dict([(k, getattr(coordvar, k)) for k in coordvar.ncattrs()])
        outf.createVariable(coordkey, coordvar.dtype.char, coordvar.dimensions, values = coordvar[:], **propd)
    return outf


def pncdump(f, outpath, format):
    from netCDF4 import Dataset
    outf = Dataset(outpath, 'w', format = format)
    for dk, dv in f.dimensions.iteritems():
        if dv.isunlimited():
            outf.createDimension(dk, None)
        else:
            outf.createDimension(dk, len(dv))
    
    for pk in f.ncattrs():
        setattr(outf, pk, getattr(f, pk))
    
    for vk in f.variables.keys():
        vv = f.variables[vk]
        outv = outf.createVariable(vk, vv.dtype.char.replace('l', 'i'), vv.dimensions)
        try:
            outv[:] = vv[:]
        except:
            # I think the layer is wrong
            outv[:] = vv[:,:outv.shape[1]]
        for pk in vv.ncattrs():
            setattr(outv, pk, getattr(vv, pk))
    
    outf.close()

def run():
    from optparse import OptionParser
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    parser = MyParser(description = """
bpch.py provides scriptable plotting and acts as a NetCDF-like library for Binary Punch files. 

For plotting, the a binary punch file must be provided. If 

For use as a library, use "from bpch import bpch" in a python script. For more information, on the bpch reader execute help(bpch) after importing as described above.

""", epilog = '')
    parser.set_usage("Usage: python bpch.py [-dpC] [-n NETCDF4_OPTION] [-s SLICE_DEF] [-r REDUCE_DEF] [-g GROUP] [-v VARIABLE] [bpchpath1 [bpchpath2 [... [bpchpathN]]]")
    parser.add_option("-d", "--difference", dest="difference",action="store_true",default=False,
                        help="Plot (bpchpath1 - bpchpath2)")

    parser.add_option("-p", "--percent-difference", dest="percent",action="store_true",default=False,
                        help="Plot bpchpath1 / bpchpath2 or, if used with -d, (bpchpath1 - bpchpath2) / bpchpath2 * 100.")
    
    parser.add_option("-g", "--group", dest = "group", default = None,
                        help = "bpch variables are organized into groups whose names are defined in diaginfo.dat; Provide a list of group names separated by ','")

    parser.add_option("-v", "--variable", dest = "variable", default = None,
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")
    
    parser.add_option("-n", "--netcdf", dest = "netcdf", default = 'NETCDF4_CLASSIC',
                        help = "NetCDF output version (options=NETCDF3_CLASSIC,NETCDF4_CLASSIC,NETCDF4; default=NETCDF4_CLASSIC).")
    
    parser.add_option("", "--vertgrid", dest = "vertgrid", default = 'GEOS-5-REDUCED',
                        help = "Vertical grid used in bpch file (options=GEOS-5-REDUCED,GEOS-5-NATIVE,MERRA-REDUCED,MERRA-NATIVE,GEOS-4-REDUCED,GEOS-4-NATIVE; default=GEOS-5-REDUCED).")
    
    parser.add_option("-s", "--slice", dest = "slice", type = "string", action = "append", default = [],
                        help = "bpch variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")
    
    parser.add_option("-r", "--reduce", dest = "reduce", type = "string", action = "append", default = [], help = "bpch variable dimensions can be reduced using dim,function,weight syntax (e.g., --reduce=layer,mean,weight). Weighting is not fully functional.")

    parser.add_option("-o", "--outpath", dest = "outpath", type = "string", action = "append", default = [], help = "bpch variable dimensions can be reduced using dim,function syntax")

    parser.add_option("-C", "--clobber", dest="clobber",action="store_true",default=False,
                        help="Overwrite pre-existing files (default = False)")


    parser.add_option("", "--oldplot", dest = "oldplot", action='store_true', default = False, help = "Use old plot interface")

    try:
        (options, args) = parser.parse_args()
    except:
        import sys
        if sys.argv[1] == '--oldplot':
            runold()
            exit()
    nfiles = len(args)
    if nfiles == 0:
        parser.print_help()
        exit()
    outfs = []
    npad = len(args) - len(options.outpath)
    if len(options.outpath) == 0:
        options.outpath.append('new_%s.nc')
    options.outpath.extend(npad * options.outpath[-1:])
    for fpath, npath in zip(args, options.outpath):
        bf = bpch(fpath, vertgrid = options.vertgrid)
        varkeys = None
        if not options.group is None:
            varkeys = reduce(list.__add__, [[(grp.strip(), k) for k in bf.groups[grp].variables.keys() if k not in metakeys] for grp in options.group.split(',')])
        else:
            varkeys = reduce(list.__add__, [[(gk,k) for k in g.variables.keys() if k not in metakeys] for gk, g in bf.groups.iteritems()])
        if not options.variable is None:
            varre = [re.compile(var) for var in options.variable.split(',')]
            varkeys = [gk.strip() + '_' + k for gk, k in varkeys if any([var.match(k) for var in varre])]
        else:
            varkeys = [gk.strip() + '_' + k for gk, k in varkeys]
        if varkeys is not None:
            varkeys.extend(['AREA', 'VOL'])
        f = getvarpnc(bf, varkeys)
        for opts in options.slice:
            f = slice_dim(f, opts)
        for opts in options.reduce:
            f = reduce_dim(f, opts)
        
        if '%s' in npath:
            newpath = npath % os.path.basename(fpath)
        else:
            newpath = npath
        
        if not options.clobber and os.path.exists(newpath):
            ask = raw_input('Enter Y to overwrite; S to skip; A to abort\n:').lower()
            if ask == 'a': exit()
            if ask == 'y': pncdump(f, newpath, format = options.netcdf)
        else:
            pncdump(f, newpath, format = options.netcdf)
        
        if options.difference or options.percent:
            outfs.append(f)
        del bf
        gc.collect()
            
    if options.difference or options.percent:
        f = getdiffpnc(*outfs, percent = options.percent)
        pncdump(f, 'diff_' + os.path.basename(args[0]) + '_' + os.path.basename(args[1]))


def runold():
    from optparse import OptionParser
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    parser = MyParser(description = """
bpch.py provides scriptable plotting and acts as a NetCDF-like library for Binary Punch files. 

For plotting, the a binary punch file must be provided. If 

For use as a library, use "from bpch import bpch" in a python script. For more information, on the bpch reader execute help(bpch) after importing as described above.

""",
    epilog="""
If -d and -p are not provided, a variable (specified by -v)
from a group (specified by -g) from each bpchpath (1...N) 
is plotted in 2 dimensions after using -t, -l, -r, and -c to 
reduce other dimensions.

Options -n -x -t -l -r -c and --title can be specified once 
for each for each bpchpath. If only one value is provided, it 
is repeated for each bpchpath

-t -l -r -c (or --time, --layer, --row, --col) are all dimension 
options. If an integer is provided, a single slice of that 
dimensions is taken using 0-based index (i.e., the first time 
is time 0). If a function is provided, it will be applied to 
that dimension (e.g., mean, median, max, min, std...). If a non-function 
value is provided, it will be converted to a python object 
and be used to acquire a "slice." A special case is "animate" which will
produce mp4 files instead of png files. "animate" requires that you have
a working installation of ffmpeg.


**Time is processed before layer, layer is processed before row, 
row is processed before col. If -t is std and -t is std:
    plotvals = timestd(layerstd(var))
  
Examples:  
    
    1. This example produce a Lat-Lon, time average, mean 
       layer with a log color-scale.

       $ python bpch.py -g IJ-AVG-$ -v O3  -t mean -l mean --log ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timemean_layermean_rowall_colall.png


    2. This example produces a zonal mean, time average plot 
       stopping at 20km (if BOXHEIGHT available) or the 21st layer.
 
       $ python bpch.py -g IJ-AVG-$ -v O3 -t mean -c mean --ymax 20 ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colmean.png

    
    3. This example produces a 1st layer Latitude-Time Hovmoller Diagram. 
    
       $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -c mean ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowall_colmean.png


    4. This example produces a 1st layer Longitude-Time Hovmoller Diagram. 
    
       $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -r mean ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowmean_colall.png

    
    5. This example would produce two Ox figures. Both from 
       time 1 (default), but the first file from layer 1 and
       the second file from layer 2. The first figure has a 
       minimum (max) value of 20 (60) and the second has a 
       minimum (maximum) of 25 (65).

       $ python bpch.py -g IJ-AVG-$ -v O3 -n 20 -x 60 -t 0 -l 0 -n 25 -x 65 -t 0 -l 1 ctm.bpch ctm.bpch2
       Successfully created ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png
       Successfully created ctm.bpch2_O3_time0_layer1_rowall_colall.png

       Successfully created ctm.bpch1_Ox_time0_layer0.png
       Successfully created ctm.bpch2_Ox_time0_layer1.png    

    6. This example would produce one Ox difference figure with 
       a minimum of -2 and a maximum of 2.

       $ python bpch.py -d -g IJ-AVG-$ -v O3 -n -2 -x 2 -t 0 -l 0 ctm.bpch ctm.bpch2
       Successfully created ctm.bpch-ctm.bpch2-diff_O3_time0_layer0_rowall_colall.png

    7. This example would produce one O3 animation across time

       $ python bpch.py -d -g IJ-AVG-$ -v O3 -t animate -l 0 ctm.bpch
       Successfully created ctm.bpch_O3_timeanimation_layer0_rowall_colall.mp4
""")
    parser.set_usage("Usage: python bpch.py [-dpnx] [-t TIME] [-l LAYER] [-r ROW] [-c COL] [-g GROUP] [-v VARIABLE] [--title=TITLE] [bpchpath1 [bpchpath2 ... bpchpathN]]")
    parser.add_option("-d", "--difference", dest="difference",action="store_true",default=False,
                        help="Plot (bpchpath1 - bpchpath2)")

    parser.add_option("-p", "--percent-difference", dest="percent",action="store_true",default=False,
                        help="Plot bpchpath1 / bpchpath2 or, if used with -d, (bpchpath1 - bpchpath2) / bpchpath2 * 100.")
    
    parser.add_option("", "--title", dest = "title", type = "string", help = "Title for plot.", action = "append", default = [])
    parser.add_option("-n", "--min", dest = "vmin", type = "float", action = 'append', default = [], help = "Minimum for the color-scale")
    parser.add_option("-x", "--max", dest = "vmax", type = "float", action = 'append', default = [], help = "Maximum for the color-scale")
    parser.add_option("", "--xmin", dest = "xmin", type = "float", action = 'append', default = [], help = "Minimum for the x-axis")
    parser.add_option("", "--xmax", dest = "xmax", type = "float", action = 'append', default = [], help = "Maximum for the x-axis")
    parser.add_option("", "--ymin", dest = "ymin", type = "float", action = 'append', default = [], help = "Minimum for the y-axis")
    parser.add_option("", "--ymax", dest = "ymax", type = "float", action = 'append', default = [], help = "Maximum for the y-axis")
    parser.add_option("", "--backend", dest = "backend", type = "string", default = 'Agg', help = "Set the backend for use. See matplotlib for more details.")
    parser.add_option("", "--log", dest = "log", action = "store_true", default = False, help = "Put data on a log scale.")
    parser.add_option("-g", "--group", dest = "group", default = None,
                        help = "bpch variables are organized into groups whose names are defined in diaginfo.dat; for a list simply do not provide the group and you will be prompted.")

    parser.add_option("-v", "--variable", dest = "variable", default = None,
                        help = "bpch variables have names defined in tracerinfo.dat; for a list simply do not provide the variable and you will be prompted.")
    
    parser.add_option("-t", "--time", dest = "time", type = "string", action = "append", default = [],
                        help = "bpch variables have time intervals and plotting currently only supports single time slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("-l", "--layer", dest = "layer", type = "string", action = "append", default = [],
                        help = "bpch variables have eta-sigma pressure layers and plotting currently only supports single layer slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("-r", "--row", dest = "row", type = "string", action = "append", default = [],
                        help = "bpch variables has rows for latitudes and plotting currently only supports single row slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("-c", "--column", dest = "column", type = "string", action = "append", default = [],
                        help = "bpch variables has rows for latitudes and plotting currently only supports single row slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("", "--format", dest="format", type="string", default="png", help = "Sets the output format (png, pdf, etc.)")

    parser.add_option("", "--mhchem", dest="mhchem", action = "store_true", default = False, help = "Use the mhchem latex options for typesetting.")

    parser.add_option("", "--latex", dest="latex", action = "store_true", default = False, help = "Use the latex options for typesetting.")

    parser.add_option("", "--dpi", dest="dpi", type="int", default=300, help = "Sets the dots per inch for output figures.")

    parser.add_option("", "--cmap", dest="cmap", type="string", default='jet', help = "Sets cmap property of tile plots. See matplotlib for more details. Good options (jet, bwr)")

    parser.add_option("", "--frames-per-second", dest="fps", type="float", default=0.5, help = "Used in combination with animate on time, layer, row, or column to regulate the number of frames per second.")

    parser.add_option("", "--oldplot", dest = "oldplot", action='store_true', default = False, help = "Use old plot interface")


    (options, args) = parser.parse_args()
    use(options.backend)
    from pylab import rcParams
    rcParams['text.usetex'] = options.latex | options.mhchem
    nfiles = len(args)
    if nfiles == 0:
        parser.print_help()
        exit()
    if options.difference or options.percent:
        fpath1 = args[0]
        fpath2 = args[1]
        f, group_key, var_key, var1 = getvar(fpath1, group_key = options.group, var_key = options.variable)
        del f
        gc.collect()
        f2, group_key2, var_key2, var2 = getvar(fpath2, group_key = options.group, var_key = options.variable)
        var_key = var1.long_name.strip() + ' - ' + var2.long_name.strip()
        if options.difference:
            var = var1[:] - var2[:]
        else:
            var = var1[:]
        if options.percent:
            var = var / var2[:] * 100
            if options.difference:
                var.units = 'pctdiff'
            else:
                var.units = 'pct'
        
        if options.difference and options.percent:
            keyappend = 'pctdiff'
        elif options.difference:
            keyappend = 'diff'
        elif options.percent:
            keyappend = 'pct'
        
        plotfvars = [('%s-%s-%s' % (os.path.basename(fpath1), os.path.basename(fpath2), keyappend), f2, '', var_key, var)]
    else:
        plotfvars = [(fpath,) + getvar(fpath, group_key = options.group, var_key = options.variable) for fpath in args]
    nplots = len(plotfvars)
    if all(map(lambda x: x == [], (options.time, options.layer, options.row, options.column))):
        warn("Using defaults for reducing dimensions: time=mean, layer = 0 (first)")
        options.time = ['mean']
        options.layer = ['0']
    times = pad(nplots, options.time, "time", "slice(None)")
    layers = pad(nplots, options.layer, "layer", "slice(None)")
    cols = pad(nplots, options.column, "column", "slice(None)")
    rows = pad(nplots, options.row, "row", "slice(None)")
    vmins = pad(nplots, options.vmin, "vmin", None)
    vmaxs = pad(nplots, options.vmax, "vmax", None)
    xmins = pad(nplots, options.xmin, "xmin", None)
    xmaxs = pad(nplots, options.xmax, "xmax", None)
    ymins = pad(nplots, options.ymin, "ymin", None)
    ymaxs = pad(nplots, options.ymax, "ymax", None)
    titles = pad(nplots, options.title, "titles", None)
    
    for time_str, layer_str, row_str, col_str, vmin, vmax, xmin, xmax, ymin, ymax, title_str, (fpath, f, group_key, var_key, var) in zip(times, layers, rows, cols, vmins, vmaxs, xmins, xmaxs, ymins, ymaxs, titles, plotfvars):
        try:
            fig_path = ('%s_%s_%s_time%s_layer%s_row%s_col%s.%s' % (os.path.basename(fpath), group_key, var_key, time_str, layer_str, row_str, col_str, options.format)).replace('-$', '').replace('$', '').replace(' ', '').replace('slice(None)', 'all')
            
            toplot = reduce_dimold(var[:], time_str, axis = 0)
            toplot = reduce_dimold(toplot, layer_str, axis = 1)
            toplot = reduce_dimold(toplot, row_str, axis = 2)
            toplot = reduce_dimold(toplot, col_str, axis = 3)
            if toplot.shape[2] != 1 and toplot.shape[3] != 1:
                maptype = 0
            elif toplot.shape[1] != 1 and toplot.shape[2] != 1:
                maptype = 1
            elif toplot.shape[1] != 1 and toplot.shape[3] != 1:
                maptype = 2
            elif toplot.shape[0] != 1 and toplot.shape[1] != 1:
                maptype = 3
            elif toplot.shape[0] != 1 and toplot.shape[2] != 1:
                maptype = 4
            elif toplot.shape[0] != 1 and toplot.shape[3] != 1:
                maptype = 5
            else:
                print toplot.shape
            
            anim_idx = None
            if 'animate' in (time_str, layer_str, row_str, col_str):
                anim_idx = [time_str, layer_str, row_str, col_str].index('animate')
            elif toplot.squeeze().ndim > 2:
                raise UserWarning("Cannot create plot with more than 2 dimensions;\n\tUse the -t, -l, -r, -c options to reduce dimensions.")
            else:
                toplot = toplot.squeeze()
            var_label = var_key
            if options.mhchem:
                var_label = r'\ce{%s}' % var_label
                
            if title_str is None:
                title_str = ('%s (Time: %s, Layer: %s, Row: %s, Col: %s)' % (var_label, time_str, layer_str, row_str, col_str)).replace('$', r'\$').replace('_', ' ').replace('slice(None)', 'all')
            if anim_idx is not None:
                from pylab import figure, draw
                from matplotlib.animation import FuncAnimation
                fig_path = fig_path.replace('.png', '.mp4')
                fig = tileplot(f, toplot[(slice(None),) * anim_idx + (0,)].squeeze(), maptype = maptype, vmin = vmin, vmax = vmax, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, title = title_str, log = options.log, cmap = options.cmap)
                def getframe(i):
                    fig.axes[0].collections[-1].set_array(toplot[(slice(None),) * anim_idx + (i,)].ravel())
                    draw()
                try:
                    animo = FuncAnimation(fig, getframe, frames = np.arange(toplot.shape[anim_idx]), interval = 1, blit = True)
                    animo.save(fig_path, codec = 'libmpeg', fps = .25)
                except IOError, e:
                    print str(e)
                    print 
                    print "-----------------------------"
                    print "Do you have ffmpeg installed?"
                    if raw_input("Should I continue creating %s files for each frame (y/n)?\n:" % options.format).lower() == 'y':
                        for i in np.arange(toplot.shape[anim_idx]):
                            getframe(i)
                            fig.savefig(fig_path.replace('.mp4', '.%d.%s' % (i, options.format)), dpi = options.dpi)
            else:
                fig = tileplot(f, toplot, maptype = maptype, vmin = vmin, vmax = vmax, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, title = title_str, log = options.log, cmap = options.cmap)
                fig.savefig(fig_path, dpi = options.dpi)
                
            f.close()
            print("Successfully created %s" % fig_path)
        except IOError, e:
            print("Unable to produce test figure (maybe you don't have matplotlib or basemap); " + str(e))
    

if __name__ == '__main__':
    run()
