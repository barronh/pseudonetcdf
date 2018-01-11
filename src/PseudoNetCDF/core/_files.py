import unittest
from PseudoNetCDF._getreader import registerreader
from PseudoNetCDF.netcdf import NetCDFFile
from collections import OrderedDict
from ._dimensions import PseudoNetCDFDimension
from ._variables import PseudoNetCDFVariable, PseudoNetCDFMaskedVariable
from warnings import warn
import numpy as np
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

class PseudoNetCDFType(type):
    """
    Create a PseudoNetCDFType meta-class
    """
    def __init__(cls, name, bases, clsdict):
        pieces = str(cls).split('\'')[1].split('.')
        longname = '.'.join([p for p in pieces[1:-1]  if '_' != p[0] and p not in ('core',)] + [pieces[-1]])
        if len(cls.mro()) > 2:
            if name not in ('PseudoNetCDFFile', 'PseudoNetCDFFileMemmap', 'WrapPnc'):
                shortl = registerreader(name, cls)
                longl = registerreader(longname, cls)
                if not (shortl or longl):
                    warn('Not registered either as ' + name + ' or ' + longname)
        super(PseudoNetCDFType, cls).__init__(name, bases, clsdict)

PseudoNetCDFSelfReg = PseudoNetCDFType('pnc', (object,), dict(__doc__ = 'Test'))

class PseudoNetCDFFile(PseudoNetCDFSelfReg):
    """
    PseudoNetCDFFile provides an interface and standard set of
    methods that a file should present to act like a netCDF file
    using the Scientific.IO.NetCDF.NetCDFFile interface.
    """
    def getMap(self, maptype = 'basemap_auto', **kwds):
        """
        maptype - choices: 'basemap', 'basemap_auto', 'cartopy' (not yet)
                  basemap - attempts to open a basemap with only the supplied kwds
                  basemap_auto - automatically adds llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat
                                 based on longitude_bounds
        """
        if maptype.startswith('basemap')
            from PseudoNetCDF.coordutil import basemap_from_proj4
            if maptype.endswith('_auto'):
                # Get edges for bounding
                lonb = f.variables['longitude_bounds']
                latb = f.variables['latitude_bounds']

                edges = dict(llcrnrlon = lonb[0, 0, 0], llcrnrlat = latb[0, 0, 0],
                             urcrnrlon = lonb[-1, -1, 2], urcrnrlat = latb[-1, -1, 2])
                kwds.update(edges)
            return basemap_from_proj4(self.getproj(withgrid = True, projformat = 'proj4'), **kwds)
        elif maptype == 'cartopy':
            raise ValueError('cartopy is not yet implemented')
        else:
            raise ValueError('maptype must be basemap, basemap_auto, or cartopy')
    
    def getproj(self, withgrid = False, projformat = 'pyproj'):
        """
        withgrid - use grid units instead of meters
        projformat - 'pyproj' (default), 'proj4' or 'wkt' allows function to return
                     a pyproj projection object or a string in the format of proj4 or WKT
        """
        if projformat == 'pyproj':
            from PseudoNetCDF.coordutil import getproj
            return getproj(self, withgrid = withgrid)
        elif projformat == 'proj4':
            from PseudoNetCDF.coordutil import getproj4
            return getproj4(self, withgrid = withgrid)
        elif projformat == 'wkt':
            from PseudoNetCDF.coordutil import getprojwkt
            return getprojwkt(self, withgrid = withgrid)
        else:
            raise ValueError('projformat must be pyproj, proj4 or wkt')

    def ll2xy(self, lon, lat):
        """
        Converts lon/lat to x distances (no false easting/northing)
        """
        return self.getproj()(lon, lat)
    
    def ll2ij(self, lon, lat):
        """
        Converts lon/lat to 0-based indicies (0,M), (0,N)
        """
        import numpy as np
        p = self.getproj(withgrid = True)
        x, y = p(lon, lat)
        i = np.asarray(x).astype('i')
        j = np.asarray(y).astype('i')
        return i, j
    
    def xy2ll(self, x, y):
        """
        Converts x, y to lon, lat (no false easting/northing)
        """
        p = self.getproj()
        lon, lat = p(x, y, inverse = True)
        return lon, lat
        
    def ij2ll(self, i, j):
        """
        Converts i, j to lon, lat (no false easting/northing)
        using cell centers assuming 0-based i/j
        """
        p = self.getproj(withgrid = True)
        lon, lat = p(i + 0.5, j + 0.5, inverse = True)
        return lon, lat
    
    def _newlike(self):
        if isinstance(self, PseudoNetCDFFile):
            outt = type(self)
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        return outf
        
    def renameVariable(self, key, inplace = False, copyall = False):
        if inplace:
            outf = self
        else:
            from ._functions import getvarpnc
            if copyall: newkeys = None
            else: newkeys = [key]
            outf = getvarpnc(self, newkeys)
        outf.variables[n] = outf.variables[o]
        del outf.variables[o]
        return outf
    def renameDimension(self, key, inplace = False):
        if inplace:
            outf = self
        else:
            from ._functions import getvarpnc
            outf = getvarpnc(self, None)
         
        outf.dimensions[n] = outf.dimensions[o]
        del outf.dimensions[o]
        for k, v in outf.variables.items():
            if o in v.dimensions:
                v.dimensions = tuple([n if d == o else d for d in v.dimensions])
        
        return outf
    
    def eval(self, expr, inplace = False, copyall = False):
        import numpy as np
        # Copy file to temporary PseudoNetCDF file
        comp = compile(expr, 'none', 'exec')
        vardict = {k: v for k, v in self.variables.items()}
        for pk in self.ncattrs():
            if not pk in vardict:
                vardict[pk] = getattr(self, pk)
        vardict['np'] = np
        vardict['self'] = self
        from symtable import symtable
        symtbl = symtable(expr, '<pncexpr>', 'exec')
        symbols = symtbl.get_symbols()
        for symbol in symbols:
            key = symbol.get_name()
            if key in vardict:
                tmpvar = vardict[key]
                break
        else:
            key = 'N/A'
            tmpvar = PseudoNetCDFVariable(None, 'temp', 'f', ())
        
        if inplace:
            outf = self
        else:
            from ._functions import getvarpnc
            if copyall: newkeys = None
            else: newkeys = [key]
            outf = getvarpnc(self, newkeys)
            try: del outf.variables[key]
            except: pass
        propd = dict([(k, getattr(tmpvar, k)) for k in tmpvar.ncattrs()])
        dimt = tmpvar.dimensions
        vardict['outf'] = self

        # Assign expression to new variable.
        exec(comp, None, vardict)
        assignedkeys = [s.get_name() for s in symbols if s.is_assigned()]
        assignedkeys = [k for k in assignedkeys if k in vardict]
        for key in assignedkeys:
            val = vardict[key]
            # if the output variable has no dimensions, there is likely a problem
            # and the output should be defined.
            if isinstance(val, (PseudoNetCDFVariable,)) and val.dimensions != ():
                outf.variables[key] = val
            else:
                outf.createVariable(key, val.dtype.char, dimt, values = val, **propd)

        return outf
    
    def setncatts(self, attdict):
        for pk, pv in attdict.items():
            setattr(self, pk, pv)
    
    def getncatts(self):
        outd = OrderedDict()
        for pk in self.ncattrs():
            outd[pk] = getattr(self, pk)
        return outd
    
    def _copywith(self, props = True, dimensions = True, variables = False, data = False):
        outf = self._newlike()
        if props: outf.setncatts(self.getncatts())
        if dimensions:
            for dk, dv in self.dimensions.items():
                newdl = len(dv)
                ndv = outf.createDimension(dk, newdl) 
                ndv.setunlimited(dv.isunlimited())
        if variables: 
            for vk, vv in self.variables.items():
                nvv = outf.createVariable(vk, vv.dtype, vv.dimensions)
                for pk in self.ncattrs():
                    pv = getattr(self, pk)
                    nvv.setncattr(pk, pv)
                    if data: nvv[:] = vv[:]
        return outf
    
    def applyAlongDimensions(self, **dimfuncs):
        """
           dimfuncs - key value pairs where the key is a dimenions and the value
                      is a 1D function (func1d) or a dictionary. If the value is a dictionary
                      it must include func1d as a function and any keyword arguments
                      as additional options
             
        """
        outf = self._copywith(props = True, dimensions = True)
        for dk, ds in dimfuncs.items():
            dv = self.dimensions[dk]
            if dk in dimfuncs:
                if dk in self.variables:
                    dvar = self.variables[dk]
                else:
                    dvar = np.arange(len(dv))
                newdl = dimfuncs[dk](dvar[:]).size
            else:
                newdl = len(dv)
            ndv = outf.createDimension(dk, newdl) 
            ndv.setunlimited(dv.isunlimited())
        
        for vark, varo in self.variables.items():
             vdims = varo.dimensions
             newvals = varo[...]
             dik = list(enumerate(vdims))
             for di, dk in dik[::-1]:
                 if dk in dimfuncs:
                     opts = dict(axis = di, arr = newvals)
                     dfunc = dimfuncs[dk]
                     if isinstance(dfunc, dict):
                         opts.update(dfunc)
                     else:
                         opts['func1d'] = dfunc
                     newvals = np.apply_along_axis(dfunc, di, newvals)
             newvaro = outf.createVariable(vark, varo.dtype, vdims)
             newvaro[...] = newvals
             for pk in varo.ncattrs():
                 setattr(newvaro, pk, getattr(varo, pk))
        
        return outf
    
    def getTimes(self):
        from PseudoNetCDF.coordutil import _parse_ref_date
        from datetime import datetime, timedelta
        if 'time' in self.variables.keys():
            time = self.variables['time']
            if 'since' in time.units:
                unit, base = time.units.strip().split(' since ')
                sdate = _parse_ref_date(base)
                out = sdate + np.array([timedelta(**{unit: float(i)}) for i in time[:]])
                return out
            else:
                return time
        elif 'TFLAG' in self.variables.keys():
            dates = self.variables['TFLAG'][:][:, 0, 0]
            times = self.variables['TFLAG'][:][:, 0, 1]
            yyyys = (dates // 1000).astype('i')
            jjj = dates % 1000
            hours = times // 10000
            minutes = times % 10000 // 100
            seconds = times % 100
            days = jjj + (hours + minutes / 60. + seconds / 3600.) / 24.
            out = np.array([datetime(yyyy, 1, 1) + timedelta(days = day - 1) for yyyy, day in zip(yyyys, days)])
            return out
        elif 'tau0' in self.variables.keys():
            out = datetime(1985, 1, 1, 0) + np.array([timedelta(hours =i) for i in self.variables['tau0'][:]])
            return out
        else:
            raise ValueError('cannot understand time for file')
    
    def subsetVariables(self, varkeys, inplace = False):    
        if inplace:
            outf = self
            for varkey in list(outf.variables):
                if not varkey in varkeys:
                    del outf.variables[varkey]
        else:
            outf = self._copywith(props = True, dimensions = True)
            for varkey in varkeys:
                varo = self.variables[varkey]
                newvaro = outf.createVariable(varkey, varo.dtype, varo.dimensions)
                for pk in varo.ncattrs():
                    setattr(newvaro, pk, getattr(varo, pk))
                newvaro[:] = varo[:]
        return outf 

    def sliceDimensions(self, newdims = ('POINTS',), **dimslices):
        """
        dimslices - key value pairs where the key is a dimension and the
                    value is a valid slice object (slices, ints or iterables)
                    if iterables are provided, all iterables must be the same
                    size and shape. If the arrays are not 1D, newdims must have ndim
                    names
        newdims - names for new dimensions. When more than one iterable applies to
                 a variable slice, fancy indexing removes both dimensions and creates
                 a new one of the iterable lengths
        """
        outf = self._copywith(props = True, dimensions = True)
        isarray = {dk: not np.isscalar(dv) for dk, dv in dimslices.items()}
        anyisarray = np.sum(list(isarray.values())) > 1
        
        if anyisarray:
            arraylens = np.array([np.asarray(da).size for dk, da in dimslices.items() if isarray[dk]])
            arrayshapes = np.array([np.asarray(da).shape for dk, da in dimslices.items() if isarray[dk]])
            arraylen = arraylens[0]
            arrayshape = arrayshapes[0]
            if not (arraylens == arraylen).all() or  not (arrayshapes == arrayshape).all():
                raise ValueError('If slicing with arrays, they must all be the same size and shape')
            for dk, ia in isarray.items():
                if ia:
                    dimslices[dk] = np.asarray(dimslices[dk])
        
        for dk, ds in dimslices.items():
            #if anyisarray and isarray[dk]: continue
            dv = self.dimensions[dk]
            if dk in dimslices:
                if dk in self.variables:
                    dvar = self.variables[dk]
                else:
                    dvar = np.arange(len(dv))
                newdl = dvar[dimslices[dk]].size
            else:
                newdl = len(dv)
            ndv = outf.createDimension(dk, newdl) 
            ndv.setunlimited(dv.isunlimited())
        
        if anyisarray:
            for ni, newdim in enumerate(newdims):
                outf.createDimension(newdim, arrayshape[ni])
        
        for vark, varo in self.variables.items():
             odims = vdims = varo.dimensions
             sliceo = tuple(dimslices.get(dk, slice(None)) for dk in vdims)
             isdarray = [isarray.get(dk, False) for dk in vdims]
             needsfancy = sum(isdarray) > 1
             concatax = np.argmax(isdarray)
             if anyisarray and needsfancy:
                  odims = [dk for dk in vdims if not isarray.get(dk, False)]
                  for newdim in newdims[::-1]:
                      odims.insert(concatax, newdim)
             
             newvaro = outf.createVariable(vark, varo.dtype, odims)
             for pk in varo.ncattrs():
                 setattr(newvaro, pk, getattr(varo, pk))
             if anyisarray and needsfancy:
                 point_arrays = []
                 for ii in range(arraylen):
                      sliceoi = []
                      for si in sliceo:
                          if np.isscalar(si):
                              sliceoi.append([si])
                          elif isinstance(si, slice):
                              sliceoi.append(si)
                          else:
                              sliceoi.append(si.ravel()[ii])
                      sliceoi = tuple(sliceoi)
                      point_arrays.append(np.expand_dims(varo[sliceoi], axis = concatax))
                 newvals = np.concatenate(point_arrays, axis = concatax)
             else:
                 newvals = varo[sliceo]
             try: newvaro[:] = newvals
             except: newvaro[:] = newvals.reshape(newvaro.shape)
        
        return outf
             
    def __repr__(self):
        from PseudoNetCDF.pncdump import pncdump
        import sys
        if sys.version_info.major == 3:
            from io import StringIO
        else:
            from StringIO import StringIO
        
        out = StringIO()
        pncdump(self, header = True, outfile = out)
        out.seek(0, 0)
        return out.read()
    
    @classmethod
    def isMine(cls, *args, **kwds):
        return False
            
    def __new__(mcl, *args, **kwds):
        new = super(PseudoNetCDFFile, mcl).__new__(mcl)
        new.variables = OrderedDict()
        new.dimensions = OrderedDict()
        new._ncattrs = ()
        new._operator_exclude_vars = ()
        return new
    
    def __init__(self, *args, **properties):
        for k, v in properties.items():
            setattr(self, k, v)

    def __setattr__(self, k, v):
        if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
            if k not in self._ncattrs:
                self._ncattrs += (k, )
        object.__setattr__(self, k, v)
    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def createDimension(self, name, length):
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
        import numpy as np
        if type == 'S': type = 'c'
        if isinstance(properties.get('values', 1), np.ma.MaskedArray) or 'fill_value' in properties:
            var = self.variables[name] = PseudoNetCDFMaskedVariable(self, name, type, dimensions, **properties)
        else:
            var = self.variables[name] = PseudoNetCDFVariable(self, name, type, dimensions, **properties)
        return var

    def close(self):
        """
        Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
        """
        pass

    def ncattrs(self):
        return self._ncattrs

    def setncattr(self, k, v):
        return setattr(self, k, v)
    
    def delncattr(self, k):
        self.__delattr__(k)
    
    def __add__(self, lhs):
        from ._functions import pncbo
        return pncbo(op = '+', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __sub__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '-', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __mul__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '*', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __div__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '/', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __floordiv__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '//', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __pow__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '**', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __and__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '&', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __or__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '|', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)

    def __xor__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '^', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)

    def __mod__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '%', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __lt__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '<', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __gt__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '>', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)

    def __eq__(self, lhs):
        if isinstance(lhs, (NetCDFFile, PseudoNetCDFFile)):
            from _functions import pncbo
            return pncbo(op = ' == ', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
        else:
            return lhs.__eq__(self)

    def __le__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '<=', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    def __ge__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '>=', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)

    def __ne__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '!=', ifile1 = self, ifile2 = lhs, verbose = 0, coordkeys = self._operator_exclude_vars)
    
    
    sync = close
    flush = close

class PseudoNetCDFFileMemmap(PseudoNetCDFFile):
    """
    Provides basic PseudoNetCDFFile functionality, but
    does not require that variables be created in memmory
    """
    def createVariable(self, name, type, dimensions, map, keep = True):
        var = PseudoNetCDFVariableMemmap(self, name, type, dimensions, map)
        if keep:
            self.variables[name] = var
        return var

class PseudoNetCDFVariables(OrderedDefaultDict):
    """
    PseudoNetCDFVariables provides a special implementation
    of the default dictionary that provides efficient access 
    to variables of a PseudoNetCDFFile.  PseudoNetCDFFiles may
    have large variables that should only be loaded if accessed.
    PseudoNetCDFVariables allows a user to specify a function
    that can create variables on demand.
    """
    def __init__(self, func, keys):
        """
        func: Function that takes a key and provides a 
              PseudoNetCDFVariable
        keys: list of keys that the dictionary should
              act as if it has
        """
        super(PseudoNetCDFVariables, self).__init__()
        self.__func = func
        self.__keys = keys
    def __missing__(self, k):
        """
        If the dictionary does not have a key, check if the
        user has provided that key.  If so, call the user 
        specifie function to create the variable.
        """
        if k in self.keys():
            return self.__func(k)
        else:
            raise KeyError('missing "%s"' % (k, ))

    def addkey(self, k):
        """
        Allow the user to extend keys after the object
        has been created.
        """
        if not k in self.__keys:
            self.__keys.append(k)

    def keys(self):
        return tuple(self.__keys + [k for k in dict.keys(self) if k not in self.__keys])

    def __len__(self):
        return len([k for k  in self.keys()])
    
    def items(self):
        return [(k, self[k]) for k in self.keys()]
    
    def __contains__(self, k):
        return k in [k for k in self.keys()]

class PseudoNetCDFTest(unittest.TestCase):
    def setUp(self):
        self.tncf = PseudoNetCDFFile()
        
    def testNetCDFFileInit(self):
        from numpy import arange
        tncf = self.tncf
        self.assert_(tncf.variables == {})
        self.assert_(tncf.dimensions == {})

        tncf.createDimension('TIME', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        self.assert_(len(tncf.dimensions['TIME']) == 24)
        self.assert_(len(tncf.dimensions['LAY']) == 4)
        self.assert_(len(tncf.dimensions['ROW']) == 5)
        self.assert_(len(tncf.dimensions['COL']) == 6)
        
        tncf.fish = 2
        setattr(tncf, 'FROG-DOG', 'HAPPY')

        o3 = tncf.createVariable('O3', 'f', ('TIME', 'LAY', 'ROW', 'COL'))
        self.assert_(tncf.variables.keys() == ['O3'])
        
        o3[:] = arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)
        o3.units = 'ppbv'
        self.assert_((o3 == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        self.assert_((tncf.variables['O3'] == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        self.assert_(o3.typecode() == 'f')

        filedims = tncf.dimensions.keys()
        filedims.sort()
        vardims = list(o3.dimensions)
        vardims.sort()

        self.assert_(filedims == vardims)
        from PseudoNetCDF.pncgen import Pseudo2NetCDF
        n = Pseudo2NetCDF().convert(tncf)
        self.assertEqual(n.variables.keys(), ['O3'])
        self.assertEqual(dict([(k, len(v)) for k, v in n.dimensions.items()]), {'TIME': 24, 'LAY': 4, 'ROW': 5, 'COL': 6})
        self.assert_((n.variables['O3'][...] == tncf.variables['O3'][...]).all())
        self.assert_(n.variables['O3'].units == 'ppbv')
        self.assertEqual(n.fish, 2)
        self.assertEqual(getattr(n, 'FROG-DOG'), 'HAPPY')
        
    def testNetCDFVariables(self):
        from numpy import arange
        tncf = PseudoNetCDFFile()
        tncf.createDimension('TIME', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        const = lambda *args, **kwds: PseudoNetCDFVariable(tncf, args[0], 'f', ('TIME', 'LAY', 'ROW', 'COL'), values = arange(24 * 4 * 5 * 6).reshape((24, 4, 5, 6)))
        tncf.variables = PseudoNetCDFVariables(const, ['NO', 'O3'])
        self.assert_((tncf.variables['O3'] == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        
    def runTest(self):
        pass
    

if __name__ == '__main__':
    unittest.main()
