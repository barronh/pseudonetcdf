from __future__ import print_function, unicode_literals
__all__ = ['bpch', 'ncf2bpch']
import os
import gc
import re

# part of the default Python distribution
from collections import OrderedDict
from warnings import warn
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

# numpy is a very common installed library
import numpy as np
from numpy import ndarray, fromfile, memmap, dtype, arange, zeros, ceil, diff, concatenate, append, pi, sin

class OrderedDefaultDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

try:
    from PseudoNetCDF.sci_var import PseudoNetCDFDimension, PseudoNetCDFVariable, PseudoNetCDFFile
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
        def __new__(mcl, *args, **kwds):
            new = super(PseudoNetCDFFile, mcl).__new__(mcl)
            new.variables=OrderedDict()
            new.dimensions=OrderedDict()
            new._ncattrs = ()
            return new
    
        def __init__(self, *args, **properties):
            for k, v in properties.items():
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

            for k,v in kwds.items():
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

class defaultdictfromkey(OrderedDefaultDict):
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

class defaultdictfromthesekeys(OrderedDefaultDict):
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
        self._keys = [k for k in keys]
        OrderedDefaultDict.__init__(self, default_factory)
    
    def __iter__(self):
        """
        Just like dictionary, but iterates through pre-defined keys
        """
        for i in self._keys:
            yield i
            
    def keys(self):
        """
        Just like dictionary, but iterates through pre-defined keys
        """
        for i in self._keys:
            yield i
    
    def values(self):
        """
        Just like dictionary, but iterates through pre-defined key values
        """
        for k in self.keys():
            yield self[k]
    
    def items(self):
        """
        Just like dictionary, but iterates through pre-defined keys and their calculated values
        """
        for k in self.keys():
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
        val = OrderedDefaultDict.__setitem__(self, key, value)
        if key not in self._keys:
            self._keys.append(key)
        return val
    
    def __delitem__(self, key):
        """
        Delete value with key and remove pre-defined key
        """
        if key in self._keys:
            ki = self._keys.index(key)
            del self._keys[ki]
        return OrderedDefaultDict.__delitem__(self, key)
    
    def pop(self, key):
        """
        Pop value with key and remove pre-defined key
        """
        val = OrderedDefaultDict.pop(self, key)
        if key in self._keys:
            ki = self._keys.index(key)
            del self._keys[ki]
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
    def __contains__(self, key):
        return key in self._keys
    
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
        except AttributeError as e:
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
    def __init__(self, parent, datamap, tracerinfo, diaginfo, keys, noscale = False, nogroup = False):
        """
        parent: NetCDF-like object to serve dimensions
        datamap: array of pre-dimensioned and datatyped values dim(tstep, i, j, k)
        tracerinfo: dictionary of tracer data keyed by ordinal
        keys: list of keys to serve
        """
        self.noscale = noscale
        self.nogroup = nogroup
        self._tracer_data = tracerinfo
        self._diag_data = diaginfo
        self._memmap = datamap
        self._parent = parent
        self._special_keys = set(metakeys)
        unique_keys = set(keys).union(self._special_keys)
        self._keys = [k for k in keys + metakeys if k in unique_keys]
        if 'BXHGHT-$_BXHEIGHT' not in keys:
            ki = self._keys.index('VOL')
            del self._keys[ki]
        self._example_key = keys[0]
        
    def __missing__(self, key):
        from ._vertcoord import geos_etai_pressure, geos_etam_pressure, geos_hyam, geos_hyai, geos_hybm, geos_hybi
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
            dims = ('time', 'tnv')
            if key == 'time':
                data = data.mean(1)
                dims = ('time',)
            
            dtype = 'd'
            kwds = dict(units = 'hours since 1985-01-01 00:00:00 UTC', base_units = 'hours since 1985-01-01 00:00:00 UTC', standard_name = key, long_name = key, var_desc = key)
            if key == 'time':
                kwds['bounds'] = 'time_bounds'
        elif key == 'hyai':
            data = self._parent.Ap
            dims = ('layer_bounds', )
            dtype = data.dtype.char
            if dims[0] not in self._parent.dimensions:
                self._parent.createDimension(dims[0], data.size)
            kwds = dict(units = "hPa", long_name = "hybrid A coefficient at layer interfaces", note = "unit consistent with GEOS-Chem pressure outputs")
        elif key == 'hyam':
            data = geos_hyam[self._parent.vertgrid]
            dims = ('layer', )
            dtype = data.dtype.char
            if dims[0] not in self._parent.dimensions:
                self._parent.createDimension(dims[0], data.size)
            kwds = dict(units = "hPa", long_name = "hybrid B coefficient at layer midpoints", note = "unit consistent with GEOS-Chem pressure outputs")
        elif key == 'hybi':
            data = self._parent.Bp
            dims = ('layer_bounds',)
            dtype = data.dtype.char
            if dims[0] not in self._parent.dimensions:
                self._parent.createDimension(dims[0], data.size)
            kwds = dict(units = "1", long_name = "hybrid B coefficient at layer interfaces")
        elif key == 'hybm':
            data = geos_hybm[self._parent.vertgrid]
            dims = ('layer', )
            dtype = data.dtype.char
            if dims[0] not in self._parent.dimensions:
                self._parent.createDimension(dims[0], data.size)
            kwds = dict(units = "1", long_name = "hybrid B coefficient at layer midpoints")
        elif key[:5] == 'layer':
            nlays = len(self._parent.dimensions['layer'])
            nedges = nlays + 1
            if key[5:] in (str(nedges), '_edges', '_bounds'):
                data = np.arange(0, nedges, dtype = 'f')
            else:
                data = np.arange(1, nedges, dtype = 'f')
            
            data = data[:len(self._parent.dimensions[key])]
            dims = (key,)
            dtype = 'f'
            kwds = dict(units = 'level', base_units = 'model layer', standard_name = 'model layer', long_name = key, var_desc = key, axis = "Z")
            kwds = dict(standard_name = "hybrid_sigma_pressure",
                        long_name = "hybrid level at layer midpoints",
                        units = "level",
                        positive = "up",)
        elif key in ('etai_pressure', 'etam_pressure'):
            nlays = len(self._parent.dimensions['layer'])
            nedges = nlays + 1
            if key[:4] == 'etai':
                data = geos_etai_pressure[self._parent.vertgrid]
            else:
                data = geos_etam_pressure[self._parent.vertgrid]
            
            data = data[:nedges]
            if 'etai' in key:
                dims = ('layer_bounds',)
            else:
                dims = ('layer',)
            dtype = 'f'
            kwds = dict(units = 'hPa', base_units = 'hPa', standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate', long_name = key, var_desc = key)
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
            key = str(key)
            header = self._memmap[key]['header'][0]
            sl, sj, si = header['f14'][::-1] - 1
            group = header['f7'].strip()
            offset = self._diag_data.get(group, {}).get('offset', 0)
            tracerid = header['f8']
            ord = header['f8'] + offset
            base_units = header['f9']
            reserved = header['f12']
            scale = self._tracer_data[ord]['SCALE']
            molwt = self._tracer_data[ord]['MOLWT']
            carbon = self._tracer_data[ord]['C']
            units = self._tracer_data[ord]['UNIT']
            tmp_data = self._memmap[key]['data']
            dims = ('time', 'layer', 'latitude', 'longitude')
            if len(['layer' in dk_ for dk_ in self._parent.dimensions]) > 1:
                dims = ('time', 'layer%d' % tmp_data.dtype['f1'].shape[0], 'latitude', 'longitude')
            kwds = dict(scale = scale, kgpermole = molwt, carbon = carbon, units = units, base_units = base_units, standard_name = key, long_name = key, var_desc = key, coordinates = ' '.join(dims), grid_mapping = "crs", reserved = reserved, tracerid = tracerid, category = group)
                
            try:
                assert((tmp_data['f0'] == tmp_data['f2']).all())
            except:
                raise ValueError('Could not parse with bpch; try bpch2')
            
            if self.noscale:
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
    print(f.dimensions)
    print(var.unit)
    print(var.dimensions)
    print(var.shape)
    
    """

    def __init__(self, bpch_path, tracerinfo = None, diaginfo = None, mode = 'r', timeslice = slice(None), noscale = False, vertgrid = 'GEOS-5-REDUCED', nogroup = False):
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
        from ._vertcoord import geos_hyai, geos_hybi, geos_eta_slp
        self._ncattrs = () 
        self.noscale = noscale
        self.nogroup = nogroup
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
        self.createDimension('tnv', 2)
        tracerinfo = tracerinfo or os.path.join(os.path.dirname(bpch_path), 'tracerinfo.dat')
        if not hasattr(tracerinfo, 'seek') or not hasattr(tracerinfo, 'readlines'):
            if not os.path.exists(tracerinfo) and tracerinfo != ' ':
                tracerinfo = 'tracerinfo.dat'
            if os.path.exists(tracerinfo):
                if os.path.isdir(tracerinfo): tracerinfo = os.path.join(tracerinfo, 'tracerinfo.dat')
            
            tracerinfo = open(tracerinfo)
        self._tracerinfofile = tracerinfo
        if hasattr(tracerinfo, 'readlines') and hasattr(tracerinfo, 'seek'):
            tracer_data = dict([(int(l[52:61].strip()), dict(NAME = l[:8].strip(), FULLNAME = l[9:39].strip(), MOLWT = float(l[39:49]), C = int(l[49:52]), TRACER = int(l[52:61]), SCALE = float(l[61:71]), UNIT = l[72:].strip())) for l in tracerinfo.readlines() if l[0] not in ('#', ' ')])
            tracer_names = dict([(k, v['NAME']) for k, v in tracer_data.items()])
        else:
            warn('Reading file without tracerinfo.dat means that names and scaling are unknown')
            tracer_data = OrderedDefaultDict(lambda: dict(SCALE = 1., C = 1., MOLWT = 1., UNIT = 'unknown', FULLNAME = 'unknown', NAME = 'unknown'))
            tracer_names = defaultdictfromkey(lambda key: key)
        
        diaginfo = diaginfo or os.path.join(os.path.dirname(bpch_path), 'diaginfo.dat')
        if not hasattr(diaginfo, 'readlines') or not hasattr(diaginfo, 'seek'):
            if not os.path.exists(diaginfo):
                diaginfo = 'diaginfo.dat'
            if os.path.exists(diaginfo):
                if os.path.isdir(diaginfo): diaginfo = os.path.join(diaginfo, 'diaginfo.dat')
            diaginfo = open(diaginfo, 'rt')

        self._diaginfofile = diaginfo

        if hasattr(diaginfo, 'read') and hasattr(diaginfo, 'seek'):
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
        self._groups = OrderedDefaultDict(set)
        while first_header is None or \
              offset < file_size:
            header = memmap(bpch_path, offset = offset, shape = (1,), dtype = _datablock_header_type, mode = mode)[0]
            
            group = header[7].decode().strip()
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
                        tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 0, UNIT = unit, MOLWT = 1., FULLNAME = 'unknown', NAME = 'unknown')
                    else:
                        tracername = tracer_names[tracer_number]
                        tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = tracer_data[tracer_number]['C'], UNIT = unit, MOLWT = 1.)

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
                    if self.nogroup:
                        keys.append('%s' % (tracername,))
                    else:
                        keys.append('%s_%s' % (group, tracername))
                
                break
            dim = header[13][::-1]
            start = header[14][::-1]
            data_type = dtype('>i4, %s>f4, >i4' % str(tuple(dim[:])))
            assert(data_type.itemsize == header[-2])
            data_types.append(data_type)
            if self.nogroup:
                keys.append('%s' % (tracername,))
            else:
                keys.append('%s_%s' % (group, tracername))

        # Repeating tracer indicates end of timestep
        keys = [str(k) for k in keys]
        time_type = dtype([(k, dtype([(str('header'), _datablock_header_type), (str('data'), d)])) for k, d in zip(keys, data_types)])
        field_shapes = set([v[0].fields['data'][0].fields['f1'][0].shape for k, v in time_type.fields.items()])
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

                outfile = open(outpath, 'w')
                infile = open(bpch_path, 'r')
                hdr = infile.read(hdrsize)
                outfile.write(hdr)
                for t in all_times:
                    if t in times:
                        outfile.write(infile.read(time_type.itemsize))
                    else:
                        infile.seek(time_type.itemsize, 1)
                outfile.close()
                #print(mint, maxt, nt, nt * time_type.itemsize)
                #cmd = 'dd if=%s ibs=%d skip=1 obs=%d | dd of=%s bs=%d skip=%d count=%d' % (bpch_path, hdrsize, time_type.itemsize, outpath, time_type.itemsize, mint, nt)
                #print(cmd)
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
        layerkeys = ['layer_bounds'] + ['layer%d' % l for l in layerns]
        keys.extend(layerkeys)
        keys.extend(['hyai', 'hyam', 'hybi', 'hybm', 'etai_pressure', 'etam_pressure'])
        
        # Ap [hPa]
        self.Ap = geos_hyai[self.vertgrid]

        # Bp [unitless]
        self.Bp = geos_hybi[self.vertgrid]

        self.createDimension('layer', self.Ap.size - 1)
        self.createDimension('layer_bounds', self.Ap.size)

        self.variables = _tracer_lookup(parent = self, datamap = datamap, tracerinfo = tracer_data, diaginfo = diag_data, keys = keys, noscale = self.noscale, nogroup = self.nogroup)
        del datamap
        tdim = self.createDimension('time', self.variables['tau0'].shape[0])
        tdim.setunlimited(True)
        self.groups = dict([(k, _diag_group(self, k, v)) for k, v in self._groups.items()])
        for grp in self.groups.values():
            for dk, d in self.dimensions.items():
                dmn = grp.createDimension(dk, len(d))
                dmn.setunlimited(d.isunlimited())
        self.Conventions = "CF-1.6"


    def __repr__(self):
        return PseudoNetCDFFile.__repr__(self) + str(self.variables)

def ncf2bpch(ncffile, outpath, verbose = False):
    outfile = open(outpath, 'wb')
    _general_header_type = np.dtype(dict(names = ['SPAD1', 'ftype', 'EPAD1', 'SPAD2', 'toptitle', 'EPAD2'], formats = '>i4, S40, >i4, >i4, S80, >i4'.split()))
    _datablock_header_type = np.dtype(dict(names = ['SPAD1', 'modelname', 'modelres', 'halfpolar', 'center180', 'EPAD1', 'SPAD2', 'category', 'tracerid', 'unit', 'tau0', 'tau1', 'reserved', 'dim', 'skip', 'EPAD2'], formats = '>i4, S20, 2>f4, >i4, >i4, >i4, >i4, S40, >i4, S40, >f8, >f8, S40, 6>i4, >i4, >i4'.split(', ')))
    varkeys = [k for k, v in ncffile.variables.items() if hasattr(v, 'tracerid') and k not in ['tau0', 'tau1', 'time', 'time_bounds', 'latitude', 'longitude', 'latitude_bounds', 'longitude_bounds', 'AREA', 'crs', 'layer', 'layer_edges', 'layer_bounds', 'hyam', 'hyai', 'etai_pressure', 'etam_pressure'] + ['layer%d' % i for i in range(200)]]
    var_types = []
    for varkey in varkeys:
        var = ncffile.variables[varkey]
        data_type = '%s>f' % str(tuple(var[0].shape))
        var_types.append(np.dtype(dict(names = ['header', 'SPAD1', 'data', 'EPAD1'], formats = [_datablock_header_type, '>i', data_type, '>i'])))
    general_header = np.zeros((1,), dtype = _general_header_type)
    general_header['SPAD1'] = general_header['EPAD1'] = 40
    general_header['SPAD2'] = general_header['EPAD2'] = 80
    general_header['ftype'] = ncffile.ftype.ljust(80)
    general_header['toptitle'] = ncffile.toptitle.ljust(80)
    general_header.tofile(outfile)
    time_data = zeros((1,), dtype = np.dtype(dict(names = varkeys, formats = var_types)))
    
    for varkey in varkeys:
        for attrk in _datablock_header_type.names[1:5]:
            time_data[varkey]['header'][attrk] = getattr(ncffile, attrk)
        
        time_data[varkey]['header']['SPAD1'] = time_data[varkey]['header']['EPAD1'] = 36
        time_data[varkey]['header']['SPAD2'] = time_data[varkey]['header']['EPAD2'] = 168
    
    for ti, (tau0, tau1) in enumerate(zip(ncffile.variables['tau0'], ncffile.variables['tau0'])):
        for varkey in varkeys:
            var = ncffile.variables[varkey]
            if not hasattr(var, 'tracerid'): continue
            if verbose: print(ti, tau0, tau1, varkey)
            vals = var[ti]
            header = time_data[varkey]['header']
            data = time_data[varkey]['data']
            header['tau0'] = tau0
            header['tau1'] = tau1
            header['reserved'] = getattr(var, 'reserved', ' ').ljust(40)
            header['tracerid'] = var.tracerid
            header['category'] = var.category.ljust(40)
            header['unit'] = var.base_units
            header['dim'] = list(vals.shape[::-1]) + [getattr(var, k_) + 1 for k_ in 'STARTI STARTJ STARTK'.split()]
            time_data[varkey]['SPAD1'] = time_data[varkey]['EPAD1'] = np.prod(vals.shape) * 4
            header['skip'] = time_data[varkey]['SPAD1'] + 8
            if ncffile.noscale:
                data[:] = vals
            else:
                data[:] = vals / var.scale
        time_data.tofile(outfile)
    
    
    outfile.flush()
    outdir = os.path.dirname(outpath)
    tracerpath = os.path.join(outdir, 'tracerinfo.dat')
    diagpath = os.path.join(outdir, 'diaginfo.dat')
    if hasattr(ncffile, '_tracerinfofile'):
        if not os.path.exists(tracerpath):
            ncffile._tracerinfofile.seek(0, 0)
            ncffile._tracerinfofile.seek(0, 0)
            outtrace = open(tracerpath, 'w')
            outtrace.write(ncffile._tracerinfofile.read())
            outtrace.flush()
    if hasattr(ncffile, '_diaginfofile'):
        if not os.path.exists(diagpath):
            ncffile._diaginfofile.seek(0, 0)
            outdiag = open(diagpath, 'w')
            outdiag.write(ncffile._diaginfofile.read())
            outdiag.flush()
    return outfile

from PseudoNetCDF._getwriter import registerwriter
registerwriter('bpch', ncf2bpch)

import unittest
class TestMemmaps(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF.testcase import geoschemfiles_paths
        self.bpchpath=geoschemfiles_paths['bpch']
    
    def testNCF2BPCH(self):
        bpchfile=bpch(self.bpchpath, noscale = True)
        outpath = self.bpchpath + '.check'
        from PseudoNetCDF.pncgen import pncgen
        pncgen(bpchfile,outpath, inmode = 'r', outmode = 'w', format = 'bpch', verbose = False)
        orig = open(self.bpchpath, 'rb').read()
        new = open(outpath, 'rb').read()
        assert(orig == new)
        os.remove(outpath)
        from PseudoNetCDF.sci_var import reduce_dim, slice_dim
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        ALD2_check = np.array([1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02, 1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02, 2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02, 3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02, 3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02, 1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02, 1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02, 2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02, 3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02, 3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02, 1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02, 1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02, 2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02, 3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02, 3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02]).reshape(ALD2.shape)
        bpchfile = slice_dim(reduce_dim(bpchfile, 'layer,mean'), 'time,0')
        pncgen(bpchfile,outpath, inmode = 'r', outmode = 'w', format = 'bpch', verbose = False)
        ALD2_check = ALD2_check[0].mean(0)[None, None]
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        np.testing.assert_allclose(ALD2, ALD2_check * 1e-9)
        bpchfile = bpch(outpath)
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        np.testing.assert_allclose(ALD2, ALD2_check)
        
         
    def testBPCH(self):
        bpchfile=bpch(self.bpchpath)
        ALD2=bpchfile.variables['IJ-AVG-$_ALD2']
        ALD2g=bpchfile.groups['IJ-AVG-$'].variables['ALD2']
        ALD2_check = np.array([1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02, 1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02, 2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02, 3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02, 3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02, 1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02, 1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02, 2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02, 3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02, 3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02, 1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02, 1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02, 2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02, 3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02, 3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02]).reshape(ALD2.shape)
        np.testing.assert_allclose(ALD2, ALD2_check)
        np.testing.assert_allclose(ALD2g, ALD2_check)
        np.testing.assert_allclose(bpchfile.variables['hyai'], np.array([0.0, 0.04804826, 6.593752, 13.1348, 19.61311, 26.09201, 32.57081, 38.98201, 45.33901, 51.69611, 58.05321, 64.36264, 70.62198, 78.83422, 89.09992, 99.36521, 109.1817, 118.9586, 128.6959, 142.91, 156.26, 169.609, 181.619, 193.097, 203.259, 212.15, 218.776, 223.898, 224.363, 216.865, 201.192, 176.93, 150.393, 127.837, 108.663, 92.36572, 78.51231, 56.38791, 40.17541, 28.36781, 19.7916, 9.292942, 4.076571, 1.65079, 0.6167791, 0.211349, 0.06600001, 0.01], 'f'))
        np.testing.assert_allclose(bpchfile.variables['hybi'], np.array([1.0, 0.984952, 0.963406, 0.941865, 0.920387, 0.898908, 0.877429, 0.856018, 0.8346609, 0.8133039, 0.7919469, 0.7706375, 0.7493782, 0.721166, 0.6858999, 0.6506349, 0.6158184, 0.5810415, 0.5463042, 0.4945902, 0.4437402, 0.3928911, 0.3433811, 0.2944031, 0.2467411, 0.2003501, 0.1562241, 0.1136021, 0.06372006, 0.02801004, 0.006960025, 8.175413e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'f'))
        np.testing.assert_allclose(bpchfile.variables['etai_pressure'], np.array([ 1.01325000e+03, 9.98050662e+02, 9.82764882e+02, 9.67479511e+02, 9.52195238e+02, 9.36910541e+02, 9.21625744e+02, 9.06342248e+02, 8.91059167e+02, 8.75776287e+02, 8.60493406e+02, 8.45211087e+02, 8.29929441e+02, 8.09555669e+02, 7.84087994e+02, 7.58621022e+02, 7.33159694e+02, 7.07698900e+02, 6.82238631e+02, 6.44053520e+02, 6.05879758e+02, 5.67705907e+02, 5.29549900e+02, 4.91400941e+02, 4.53269420e+02, 4.15154739e+02, 3.77070069e+02, 3.39005328e+02, 2.88927351e+02, 2.45246173e+02, 2.08244245e+02, 1.76930008e+02, 1.50393000e+02, 1.27837000e+02, 1.08663000e+02, 9.23657200e+01, 7.85123100e+01, 5.63879100e+01, 4.01754100e+01, 2.83678100e+01, 1.97916000e+01, 9.29294200e+00, 4.07657100e+00, 1.65079000e+00, 6.16779100e-01, 2.11349000e-01, 6.60000100e-02, 1.00000000e-02]))

    def runTest(self):
        pass

if __name__ == '__main__':
    run()
