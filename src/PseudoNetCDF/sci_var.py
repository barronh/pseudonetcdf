__doc__ = r"""
Scientific data has dimensions that have physical meaning and values
only have meaning in the context of their units.  This module implements
numpy arrays that are aware of their dimensions trying to vaguely adhere
to the Common Data Model from Unitdata at UCAR.

Each variable has as a property its dimensions names (dimensions).  Further,
each dimension name exists as a property and contains a one dimensional array
of values associated with that dimension.

For the purposes of ease of use, the standard properties of netCDF files
are attached and the arrays implement the Scientific.IO.NetCDF.NetCDFVariable
interfaces.
"""

__all__ = ['PseudoNetCDFFile', 'PseudoNetCDFDimension', 'PseudoNetCDFVariableConvertUnit', 'PseudoNetCDFFileMemmap', 'PseudoNetCDFVariable', 'PseudoNetCDFMaskedVariable', 'PseudoIOAPIVariable', 'PseudoNetCDFVariables', 'Pseudo2NetCDF', 'reduce_dim', 'slice_dim', 'getvarpnc', 'interpvars', 'extract', 'pncbo', 'seqpncbo']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import numpy as np
from numpy import array, zeros, ndarray, isscalar, abs
from numpy.ma import MaskedArray, zeros as mazeros
from collections import defaultdict, OrderedDict

from warnings import warn
    
from numpy import arange
from tempfile import NamedTemporaryFile as tnf
from types import MethodType
from units import convert
import re
import unittest
from os.path import isdir
import sys

from _getreader import registerreader

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
        
def PseudoNetCDFVariableConvertUnit(var,outunit):
    """
    Convert the unit of var and update the 
    associated IOAPI metadata
    """
    do = PseudoNetCDFFile()
    shape=var.shape
    for i,d in enumerate(var.dimensions):
        do.createDimension(d, shape[i])
    outvar=PseudoNetCDFVariable(do,var.long_name.strip(),var.typecode(),var.dimensions,values=convert(var,var.units,outunit))
    for k in var.ncattrs():
        v = getattr(var, k)
        setattr(outvar,k,v)
    outvar.units=outunit
    return outvar
    
class classreg(type):
    def __init__(cls, name, bases, clsdict):
        if len(cls.mro()) > 2:
            if name != 'PseudoNetCDFFileMemmap':
                registerreader(name, cls)
        super(classreg, cls).__init__(name, bases, clsdict)

class PseudoNetCDFFile(object):
    """
    PseudoNetCDFFile provides an interface and standard set of
    methods that a file should present to act like a netCDF file
    using the Scientific.IO.NetCDF.NetCDFFile interface.
    """
    __metaclass__ = classreg
    @classmethod
    def isMine(cls, *args, **kwds):
        try:
            cls(*args, **kwds)
            return True
        except:
            return False
            
    def __new__(cls, *args, **kwds):
        new = object.__new__(cls, *args, **kwds)
        new.variables = OrderedDict()
        new.dimensions = OrderedDict()
        new._ncattrs = ()
        return new
    
    def __init__(self, *args, **properties):
        for k, v in properties.iteritems():
            setattr(self, k, v)

    def __setattr__(self, k, v):
        if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
            self._ncattrs += (k,)
        object.__setattr__(self, k, v)
    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def createDimension(self,name,length):
        """
        name - string name for dimension
        length - maximum length of dimension
        """
        dim = self.dimensions[name]=PseudoNetCDFDimension(self, name, length)
        return dim

    def createVariable(self, name, type, dimensions, **properties):
        """
        name - string
        type - numpy dtype code (e.g., 'f', 'i', 'd')
        dimensions - tuple of dimension keys that can be
                     found in objects' dimensions dictionary
        """
        if isinstance(properties.get('values', 1), np.ma.MaskedArray) or 'fill_value' in properties:
            var = self.variables[name] = PseudoNetCDFMaskedVariable(self,name,type,dimensions, **properties)
        else:
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

class PseudoNetCDFFileMemmap(PseudoNetCDFFile):
    """
    Provides basic PseudoNetCDFFile functionality, but
    does not require that variables be created in memmory
    """
    def createVariable(self,name,type,dimensions,map,keep=True):
        var=PseudoNetCDFVariableMemmap(self,name,type,dimensions,map)
        if keep:
            self.variables[name]=var
        return var

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
    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)
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
        if not hasattr(self, '_ncattrs'):
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

class PseudoNetCDFMaskedVariable(PseudoNetCDFVariable, MaskedArray):
    def __setattr__(self, k, v):
        """
        Set attributes (aka properties) and identify user-defined attributes.
        """
        if not hasattr(self, k) and k[:1] != '_':
            self._ncattrs += (k,)
        MaskedArray.__setattr__(self, k, v)

    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def ncattrs(self):
        """
        Returns a tuple of attributes that have been user defined
        """
        
        return self._ncattrs

    def swapaxes(self, a1, a2):
        out = np.ma.MaskedArray.swapaxes(self, a1, a2)
        out = self.__array_wrap__(out)
        return out

    def reshape(self, shape):
        out = np.ma.MaskedArray.reshape(self, shape)
        out = self.__array_wrap__(out)
        return out
    
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
        __array_priority__ = 1000000000.
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

            result=mazeros(shape,typecode)

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
        MaskedArray.__array_finalize__(self, obj)
        if not hasattr(self, '_ncattrs'):
            self._ncattrs = ()
            if obj is None: return
            if hasattr(obj, '_ncattrs'):
                for k in obj._ncattrs:
                    setattr(self, k, getattr(obj, k))
        if not hasattr(self, 'dimensions'):
            if hasattr(obj, 'dimensions'):
                setattr(self, 'dimensions', obj.dimensions)
                self._ncattrs = self._ncattrs[:-1]

    def __array_wrap__(self, obj, context = None):
        MaskedArray.__array_finalize__(self, obj)
        if hasattr(self, '_ncattrs'):
            for k in self._ncattrs:
                setattr(obj, k, getattr(self, k))
        else:
            self._ncattrs = ()
        
        obj.dimensions = self.dimensions
        return obj
        
    def __getitem__(self, item):
        out = np.ma.MaskedArray.__getitem__(self, item)
        out.dimensions = self.dimensions
        if hasattr(self, '_ncattrs'):
            for k in self._ncattrs:
                setattr(out, k, getattr(self, k))
        else:
            self._ncattrs = ()
        return out

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

    
def PseudoIOAPIVariable(parent,name,typecode,dimensions,**kwds):
    """
    Creates a variable using the dimensions as defined in
    the parent object
    
    parent: an object with a dimensions variable
    name: name for variable
    typecode: numpy style typecode
    dimensions: a typle of dimension names to be used from
                parent
    units: default = none
    long_name: default = name
    var_desc: default = name
    """

    retval = PseudoNetCDFVariable(parent, name, typecode, dimensions, **kwds)

    if not kwds.has_key('units'):
        warn('IOAPI variables must have units; %s has been initialized with "None" units')
        retval.units = 'None'
        
    if not kwds.has_key('long_name'):
        retval.long_name = name.ljust(16)

    if not kwds.has_key('var_desc'):
        retval.var_desc = name.ljust(80)

    return retval
    
class PseudoNetCDFVariables(OrderedDefaultDict):
    """
    PseudoNetCDFVariables provides a special implementation
    of the default dictionary that provides efficient access 
    to variables of a PseudoNetCDFFile.  PseudoNetCDFFiles may
    have large variables that should only be loaded if accessed.
    PseudoNetCDFVariables allows a user to specify a function
    that can create variables on demand.
    """
    def __init__(self,func,keys):
        """
        func: Function that takes a key and provides a 
              PseudoNetCDFVariable
        keys: list of keys that the dictionary should
              act as if it has
        """
        super(PseudoNetCDFVariables, self).__init__()
        self.__func=func
        self.__keys=keys
    def __missing__(self,k):
        """
        If the dictionary does not have a key, check if the
        user has provided that key.  If so, call the user 
        specifie function to create the variable.
        """
        if k in self.keys():
            return self.__func(k)
        else:
            raise KeyError, 'missing "%s"' % (k,)

    def addkey(self,k):
        """
        Allow the user to extend keys after the object
        has been created.
        """
        self.__keys.append(k)

    def keys(self):
        return tuple(set(dict.keys(self) + self.__keys))
    
    def has_key(self,k):
        return k in self.keys()
    
    def iteritems(self):
        for k in self.keys():
            yield k,self[k]
    
    def iterkeys(self):
        for k in self.keys():
            yield k

def get_ncf_object(path_or_object, mode, format = 'NETCDF4_CLASSIC'):
    from os.path import exists, isfile
    from netcdf import NetCDFFile
    read_only = ('r', 'r+', 'rs', 'rs+', 'r+s')
    if isinstance(path_or_object, str):
        if exists(path_or_object):
            if isfile(path_or_object):
                ncf_object = NetCDFFile(path_or_object, mode, format = format)
            elif isdir(path_or_object):
                raise ValueError("Got directory at %s; not sure what to do" % path_or_object)
            else:
                raise ValueError("Expected file or directory at %s" % path_or_object)
        elif mode not in read_only:
            ncf_object = NetCDFFile(path_or_object, mode, format = format)
        else:
            raise IOError("Cannot open missing file for reading")
    elif isinstance(path_or_object, NetCDFFile) or isinstance(path_or_object, PseudoNetCDFFile):
        return path_or_object
    elif path_or_object is None and mode not in read_only:
        tfile=tnf(mode='w+b')
        npath=tfile.name
        ncf_object=NetCDFFile(npath,mode)
    else:
        raise ValueError("Not a path; not a netCDF file; not a PseudoNetCDF file... I don't know what to do")
    return ncf_object

def get_dimension_length(pfile, key):
    dim = pfile.dimsensions[key]
    if dim is None:
        for k in pfile.variables.keys():
            v = pfile.variables[k]
            if key in v.dimensions:
                return v[...].shape[list(v.dimensions).index(key)]
        return 0
    elif isinstance(dim, int):
        return dim
    else:
        return len(dim)

def interpvars(f, weights, dimension):
    """
    f - PseudoNetCDFFile
    weights - weights for new dimensions from old dimension dim(new, old)
    dimension - which dimensions will be reduced
    """
    outf = f
    if hasattr(outf, 'groups'):
        for grpk, grpv in outf.groups.items():
            outf.groups[grpk] = interpvars(grpv, weights, dimension)
    
    oldd = outf.dimensions[dimension]
    newd = outf.createDimension(dimension, weights.shape[0])
    newd.setunlimited(oldd.isunlimited())
    for vark, oldvar in outf.variables.iteritems():
        if dimension in oldvar.dimensions:
            dimidx = list(oldvar.dimensions).index(dimension)
            newvar = outf.createVariable(vark, oldvar.dtype.char, oldvar.dimensions)
            for ak in oldvar.ncattrs():
                setattr(newvar, ak, getattr(oldvar, ak))
            weightslice = (None,) * (dimidx) + (Ellipsis,) + (None,) * len(oldvar.dimensions[dimidx + 1:])
            varslice = (slice(None,),) * dimidx + (None,)
            newvar[:] = (weights[weightslice] * oldvar[varslice]).sum(dimidx + 1)
    return outf

    

def extract(f, lonlat, unique = False, gridded = None):
    from StringIO import StringIO
    outf = f
    if hasattr(outf, 'groups'):
        for grpk, grpv in outf.groups.items():
            outf.groups[grpk] = extract(grpv, lonlat)
    
    longitude = f.variables['longitude'][:]
    latitude = f.variables['latitude'][:]
    if gridded is None:
        gridded = ('longitude' in f.dimensions and 'latitude' in f.dimensions) or \
                  ('COL' in f.dimensions and 'ROW' in f.dimensions) or \
                  ('x' in f.dimensions and 'y' in f.dimensions)
    latlon1d = longitude.ndim == 1 and latitude.ndim == 1
    if latlon1d and gridded:
        latitude = latitude[(slice(None), None, None)]
        longitude = longitude[(None, slice(None), None)]
    else:
        latitude = latitude[Ellipsis, None]
        longitude = longitude[Ellipsis, None]
    
    lonlatdims = latitude.ndim - 1
    outf.lonlatcoords = ('/'.join(lonlat))
    lons, lats = np.genfromtxt(StringIO(outf.lonlatcoords.replace('/', '\n')), delimiter = ',').T
    londists = longitude - lons[(None,) * lonlatdims]
    latdists = latitude - lats[(None,) * lonlatdims]
    totaldists = ((latdists**2 + londists**2)**.5)
    
    if latlon1d and not gridded:
        latidxs, = lonidxs, = np.unravel_index(totaldists.reshape(-1, latdists.shape[-1]).argmin(0), totaldists.shape[:-1])
    else:
        latidxs, lonidxs = np.unravel_index(totaldists.reshape(-1, latdists.shape[-1]).argmin(0), totaldists.shape[:-1])
    
    if unique:
        tmpx = OrderedDict()
        for lon, lat, lonlatstr in zip(lonidxs, latidxs, outf.lonlatcoords.split('/')):
            if (lon, lat) not in tmpx:
                tmpx[(lon, lat)] = lonlatstr
        
        lonidxs, latidxs = np.array(tmpx.keys()).T
        outf.lonlatcoords_orig = outf.lonlatcoords
        outf.lonlatcoords = '/'.join([tmpx[k] for k in zip(lonidxs, latidxs)])
#     longitude = f.variables['longitude'][:]
#     latitude = f.variables['latitude'][:]
#     latidxs = []
#     lonidxs = []
#     for lon, lat in zip(lons, lats):
#         londist = lon - longitude
#         latdist = lat - latitude
#         totaldist = (latdist[:, None]**2 + londist[None, :]**2)**.5
#         latidxa, lonidxa = np.where(totaldist.min() == totaldist)
#         if len(latidxa) > 1:
#             warn("Selecting first of equidistant points")
#             
#         latidx = latidxa[0]
#         lonidx = lonidxa[0]
#         #lonidx = abs(londist).argmin()
#         #latidx = abs(latdist).argmin()
#         #import pdb; pdb.set_trace()
#         latidxs.append(latidx)
#         lonidxs.append(lonidx)
#     latidxs = array(latidxs)
#     lonidxs = array(lonidxs)
    for k, v in outf.variables.items():
        try:
            coords = v.coordinates.split()
        except:
            coords = v.dimensions
        dims = v.dimensions
        f.createDimension('points', len(latidxs))
        if 'longitude' in coords or 'latitude' in coords:
            try:
                del outf.variables[k]
            except:
                pass
            newdims = []
            for d in dims:
                if d not in ('longitude', 'latitude'):
                    newdims.append(d)
                else:
                    if 'points' not in newdims:
                        newdims.append('points')
                        
            
            newdims = tuple(newdims)
            newslice = tuple([{'latitude': latidxs, 'longitude': lonidxs, 'points': latidxs, 'PERIM': latidxs}.get(d, slice(None)) for d in dims]) 
            
            nv = outf.createVariable(k, v.dtype.char, newdims, values = v[newslice])
            for ak in v.ncattrs():
                setattr(nv, ak, getattr(v, ak))
            setattr(nv, 'coordinates', getattr(v, 'coordinates', ' '.join(coords)))
            for di, dk in enumerate(newdims):
                if dk not in outf.dimensions:
                    outf.createDimension(dk, nv.shape[di])
    return f
    
def mask_vals(f, maskdef, metakeys = 'time layer level latitude longitude time_bounds latitude_bounds longitude_bounds ROW COL LAY TFLAG ETFLAG'.split()):
    for varkey, var in f.variables.iteritems():
        if varkey not in metakeys:
            vout = eval('np.ma.masked_%s(var[:], %s)' % tuple(maskdef.split(',')))
            f.variables[varkey] = PseudoNetCDFMaskedVariable(f, varkey, var.dtype.char, var.dimensions, values = vout, **dict([(pk, getattr(var, pk)) for pk in var.ncattrs()]))
    return f
    
def slice_dim(f, slicedef, fuzzydim = True):
    """
    variables have dimensions (e.g., time, layer, lat, lon), which can be subset using 
        slice_dim(f, 'dim,start,stop,stride')
        
    e.g., slice_dim(f, 'layer,0,47,5') would sample every fifth layer starting at 0
    """
    historydef = "slice_dim(f, %s, fuzzydim = %s); " % (slicedef, fuzzydim)
    slicedef = slicedef.split(',')
    slicedef = [slicedef[0]] + map(eval, slicedef[1:])
    if len(slicedef) == 2:
        slicedef.append(slicedef[-1] + 1)
    slicedef = (slicedef + [None,])[:4]
    dimkey, dmin, dmax, dstride = slicedef    
    unlimited = f.dimensions[dimkey].isunlimited()
    if fuzzydim:
        partial_check = [key for key in f.dimensions if dimkey == key[:len(dimkey)] and key[len(dimkey):].isdigit()]
        for dimk in partial_check:
            f = slice_dim(f, '%s,%s,%s,%s' % (dimk, dmin, dmax, dstride))
        
    for varkey in f.variables.keys():
        var = f.variables[varkey]
        if dimkey not in var.dimensions:
            continue
        axis = list(var.dimensions).index(dimkey)
        vout = var[...].swapaxes(0, axis)[dmin:dmax:dstride].swapaxes(0, axis)
        
        newlen = vout.shape[axis]
        newdim = f.createDimension(dimkey, newlen)
        newdim.setunlimited(unlimited)
        f.variables[varkey] = vout
    history = getattr(f, 'history', '')
    history += historydef
    setattr(f, 'history', history)

    return f
    
def reduce_dim(f, reducedef, fuzzydim = True, metakeys = 'time layer level latitude longitude time_bounds latitude_bounds longitude_bounds ROW COL LAY TFLAG ETFLAG'.split()):
    """
    variable dimensions can be reduced using
    
    reduce_dim(file 'dim,function,weight')
    
    e.g., reduce_dim(layer,mean,weight).
    
    Weighting is not fully functional.
    """
    metakeys = [k for k in metakeys if k in f.variables.keys()]
    historydef = "reduce_dim(f, %s, fuzzydim = %s, metakeys = %s); " % (reducedef, fuzzydim, metakeys)
    import numpy as np
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
    
    unlimited = f.dimensions[dimkey].isunlimited()
    f.createDimension(dimkey, 1)
    if unlimited:
        f.dimensions[dimkey].setunlimited(True)

    for varkey in f.variables.keys():
        var = f.variables[varkey]
        if dimkey not in var.dimensions:
            continue
        
        axis = list(var.dimensions).index(dimkey)
        
        if not varkey in metakeys:
            if numweightkey is None:
                vout = getattr(var[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis)
            elif denweightkey is None:
                wvar = var * np.array(numweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)]
                vout = getattr(wvar[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis)
                vout.units = vout.units.strip() + ' * ' + numweight.units.strip()
                if hasattr(vout, 'base_units'):
                    vout.base_units = vout.base_units.strip() + ' * ' + numweight.base_units.strip()
            else:
                nwvar = var * np.array(numweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)]
                vout = getattr(nwvar[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis) / getattr(np.array(denweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]), None)], func)(axis = axis)
        else:
            if '_bounds' not in varkey and '_bnds' not in varkey:
                vout = getattr(var[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis)
            else:
                vout = getattr(var[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis)
                vout[0] = var[...].min(), var[...].max()
        f.variables[varkey] = PseudoNetCDFMaskedVariable(f, varkey, var.dtype.char, var.dimensions, values = vout)

    history = getattr(f, 'history', '')
    history += historydef
    setattr(f, 'history', history)
    return f

def pncbo(op, ifile1, ifile2):
    p2p = Pseudo2NetCDF()
    tmpfile = PseudoNetCDFFile()
    p2p.convert(ifile1, tmpfile)
    for k in tmpfile.variables.keys():
        outvar = tmpfile.variables[k]
        in1var = ifile1.variables[k]
        in2var = ifile2.variables[k]
        if outvar.ndim > 0:
            outvar[:] = eval('in1var[:] %s in2var[:]' % op)
        else:
            outvar.itemset(eval('in1var %s in2var' % op))
    return tmpfile

def seqpncbo(ops, ifiles):
    for op in ops:
        ifile1, ifile2 = ifiles[:2]
        newfile = pncbo(op, ifile1, ifile2)
        del ifiles[:2]
        ifiles.insert(0, newfile)
    return ifiles
    

def getvarpnc(f, varkeys, coordkeys = []):
    coordkeys = set(coordkeys)
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
    for propkey in f.ncattrs():
        setattr(outf, propkey, getattr(f, propkey))
    for varkey in varkeys:
        try:
            var = eval(varkey, None, f.variables)
        except:
            var = f.variables[varkey]
        for dimk, dimv in zip(var.dimensions, var.shape):
            if dimk not in outf.dimensions:
                newdimv = outf.createDimension(dimk, dimv)
                if f.dimensions[dimk].isunlimited():
                    newdimv.setunlimited(True)
                coordkeys.add(dimk)
                try:
                    tempv = f.variables[dimk]
                    if hasattr(tempv, 'bounds'):
                        coordkeys.add(tempv.bounds.strip())
                except (ValueError, KeyError, AttributeError), e:
                    pass
        for coordk in coordkeys:
            if coordk in f.dimensions and coordk not in outf.dimensions:
                newdimv = outf.createDimension(coordk, len(f.dimensions[coordk]))
                if f.dimensions[coordk].isunlimited():
                    newdimv.setunlimited(True)
    
        propd = dict([(k, getattr(var, k)) for k in var.ncattrs()])
        outf.createVariable(varkey, var.dtype.char, var.dimensions, values = var[...], **propd)
    for coordkey in coordkeys:
        if coordkey in f.variables.keys():
            coordvar = f.variables[coordkey]
            propd = dict([(k, getattr(coordvar, k)) for k in coordvar.ncattrs()])
            outf.createVariable(coordkey, coordvar.dtype.char, coordvar.dimensions, values = coordvar[...], **propd)
            for dk in coordvar.dimensions:
                if dk not in outf.dimensions:
                    dv = outf.createDimension(dk, len(f.dimensions[dk]))
                    dv.setunlimited(f.dimensions[dk].isunlimited())
    return outf



class Pseudo2NetCDF:
    """
    Pseudo2NetCDF is a base class for conversion.  Properties and methods can
    be overwritten to facilitate conversion of special PseudoNetCDFFiles.
    
    Specifically: ignore_global_properties and ignore_variable_properties lists
    can be overwritten so that class properties and methods are not written
    to a netCDF file
    """
    ignore_global_re=re.compile('^_\w*(__\w*)?')
    ignore_variable_re=re.compile('^_\w*(__\w*)?')
    ignore_global_properties=['variables','dimensions']
    ignore_variable_properties=['typecode','dimensions']
    special_properties = ['_fillvalue', '_FillValue']
    unlimited_dimensions = []
    create_variable_kwds = {}
    def __init__(self, datafirst = False, verbose = True):
        self.datafirst = datafirst
        self.verbose = verbose
    def convert(self,pfile,npath=None, inmode = 'r', outmode = 'w', format = 'NETCDF4'):
        pfile = get_ncf_object(pfile, inmode)
        nfile = get_ncf_object(npath, outmode, format = format)
        if self.verbose: print >> sys.stdout, "Adding dimensions"
        self.addDimensions(pfile,nfile)
        if self.verbose: print >> sys.stdout, "Adding globals"
        self.addGlobalProperties(pfile,nfile)
        if self.verbose: print >> sys.stdout, "Adding variables"
        self.addVariables(pfile,nfile)
        nfile.sync()
        return nfile
        
    def addDimensions(self,pfile,nfile):
        for d,v in pfile.dimensions.iteritems():
            unlim = (d in self.unlimited_dimensions or v.isunlimited())
            if not isinstance(v, (int, long)) and v is not None:
                v = len(v)
            
            if unlim:
                if isinstance(nfile, PseudoNetCDFFile):
                    nd = nfile.createDimension(d,v)
                    nd.setunlimited(True)
                else:
                    nd = nfile.createDimension(d,None)
            else:
                nd = nfile.createDimension(d,v)

        nfile.sync()
    
    def addGlobalProperties(self,pfile,nfile):
        for k in [k for k in pfile.ncattrs() if k not in self.ignore_global_properties and self.ignore_global_re.match(k)==None]:
            value=getattr(pfile,k)
            if not isinstance(value, MethodType):
                try:
                    setattr(nfile,k,value)
                except:
                    warn("Could not add %s to file" % k)

    def addVariableProperties(self,pvar,nvar):
        for a in [k for k in pvar.ncattrs() if (k not in self.ignore_variable_properties and self.ignore_variable_re.match(k)==None) or k in self.special_properties]:
            value=getattr(pvar,a)
            if not isinstance(value, MethodType):
                try:
                    setattr(nvar,a,value)
                except:
                    if 'long_name' in pvar.ncattrs():
                        warn("Could not add %s=%s to variable %s" % (a,str(value),str(pvar.long_name)))
                    else:
                        warn("Could not add %s=%s to variable" % (a,str(value)))
                        
    
    def addVariable(self,pfile,nfile,k, data = True):
        pvar=pfile.variables[k]
        try:
            typecode = pvar.typecode()
        except:
            typecode = pvar[...].dtype.char
        
        create_variable_kwds = self.create_variable_kwds.copy()
        if hasattr(pvar, 'fill_value'):
            create_variable_kwds['fill_value'] = pvar.fill_value
        nvar=nfile.createVariable(k,typecode,pvar.dimensions, **create_variable_kwds)
        self.addVariableProperties(pvar,nvar)
        if data:
            self.addVariableData(pfile, nfile, k)
        nfile.sync()
        try:
            nfile.flush()
        except:
            pass
        del pvar,nvar

    def addVariableData(self, pfile, nfile, k):
        nvar = nfile.variables[k]
        pvar = pfile.variables[k]
        if isscalar(nvar) or nvar.ndim == 0:
            nvar.assignValue(pvar)
        elif isinstance(pvar[...], MaskedArray):
            nvar[:] = pvar[...].filled()
        else:
            nvar[:] = pvar[...]
        

    def addVariables(self,pfile,nfile):
        for k in pfile.variables.keys():
            if self.verbose: print >> sys.stdout, "Defining", k
            self.addVariable(pfile,nfile,k, data = self.datafirst)
        nfile.sync()
        if not self.datafirst:
            for k in pfile.variables.keys():
                if self.verbose: print >> sys.stdout, "Populating", k
                self.addVariableData(pfile,nfile,k)
            nfile.sync()

class PseudoNetCDFTest(unittest.TestCase):
    def setUp(self):
        self.tncf=PseudoNetCDFFile()
        
    def testNetCDFFileInit(self):
        tncf=self.tncf
        self.assert_(tncf.variables=={})
        self.assert_(tncf.dimensions=={})

        tncf.createDimension('TIME',24)
        tncf.createDimension('LAY',4)
        tncf.createDimension('ROW',5)
        tncf.createDimension('COL',6)

        self.assert_(len(tncf.dimensions['TIME'])==24)
        self.assert_(len(tncf.dimensions['LAY'])==4)
        self.assert_(len(tncf.dimensions['ROW'])==5)
        self.assert_(len(tncf.dimensions['COL'])==6)
        
        tncf.fish=2
        setattr(tncf,'FROG-DOG','HAPPY')

        o3=tncf.createVariable('O3','f',('TIME','LAY','ROW','COL'))
        self.assert_(tncf.variables.keys()==['O3'])
        
        o3[:] = arange(24*4*5*6).reshape(24,4,5,6)
        o3.units='ppbv'
        self.assert_((o3==arange(24*4*5*6).reshape(24,4,5,6)).all())
        self.assert_((tncf.variables['O3']==arange(24*4*5*6).reshape(24,4,5,6)).all())
        
        self.assert_(o3.typecode()=='f')

        filedims=tncf.dimensions.keys()
        filedims.sort()
        vardims=list(o3.dimensions)
        vardims.sort()

        self.assert_(filedims==vardims)
        
        n=Pseudo2NetCDF().convert(tncf)
        self.assertEqual(n.variables.keys(),['O3'])
        self.assertEqual(dict([(k, len(v)) for k, v in n.dimensions.iteritems()]),{'TIME': 24, 'LAY': 4, 'ROW': 5, 'COL': 6})
        self.assert_((n.variables['O3'][...]==tncf.variables['O3'][...]).all())
        self.assert_(n.variables['O3'].units=='ppbv')
        self.assertEqual(n.fish,2)
        self.assertEqual(getattr(n,'FROG-DOG'),'HAPPY')
        
    def testNetCDFVariables(self):
        tncf=PseudoNetCDFFile()
        tncf.createDimension('TIME',24)
        tncf.createDimension('LAY',4)
        tncf.createDimension('ROW',5)
        tncf.createDimension('COL',6)

        const=lambda *args,**kwds: PseudoNetCDFVariable(tncf,args[0],'f',('TIME','LAY','ROW','COL'),values=arange(24*4*5*6).reshape((24,4,5,6)))
        tncf.variables=PseudoNetCDFVariables(const,['NO','O3'])
        self.assert_((tncf.variables['O3']==arange(24*4*5*6).reshape(24,4,5,6)).all())
        
        
    def runTest(self):
        pass

if __name__ == '__main__':
    unittest.main()
