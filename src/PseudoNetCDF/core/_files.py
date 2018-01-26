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
            if name not in ('PseudoNetCDFFile', 'WrapPnc'):
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
        Description
        
        Parameters
        ----------
        maptype : choices : 'basemap', 'basemap_auto', 'cartopy' (not yet)
                  basemap : attempts to open a basemap with only the supplied kwds
                  basemap_auto : automatically adds llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat
                                 based on longitude_bounds
        kwds : keywords for basemap or cartopy
        
        Returns
        -------
        map : type is basemap or cartopy axis
        """
        if maptype.startswith('basemap'):
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
        Description
        
        Parameters
        ----------
        withgrid : use grid units instead of meters
        projformat : 'pyproj' (default), 'proj4' or 'wkt' allows function to return
                     a pyproj projection object or a string in the format of proj4 or WKT
        
        Returns
        -------
        proj : string (wkt, proj4) or pyprojProj (pyproj)
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
        
        Parameters
        ----------
        lon : scalar or iterable of longitudes in decimal degrees
        lat : scalar or iterable of latitudes in decimal degrees
        
        Returns
        -------
        x, y coordinates in map projection (meters or radians)
        """
        return self.getproj()(lon, lat)
    
    def ll2ij(self, lon, lat):
        """
        Converts lon/lat to 0-based indicies (0,M), (0,N)
        
        Parameters
        ----------
        lon : scalar or iterable of longitudes in decimal degrees
        lat : scalar or iterable of latitudes in decimal degrees
        
        Returns
        -------
        i, j : indices (0-based) for variables
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
        
        Parameters
        ----------
        x : scalar or iterable of projected west-east coordinates
        y : scalar or iterable of projected south-north coordinates
        
        Returns
        -------
        lon, lat : scalars or iterables of longitudes and latitudes in decimal degrees
        """
        p = self.getproj()
        lon, lat = p(x, y, inverse = True)
        return lon, lat
        
    def ij2ll(self, i, j):
        """
        Converts i, j to lon, lat (no false easting/northing)
        using cell centers assuming 0-based i/j
        
        Parameters
        ----------
        i : scalar or iterable of indicies (0-based) for the west-east dimension
        j : scalar or iterable of indicies (0-based) for the south-north dimension
        
        Returns
        -------
        lon, lat : scalars or iterables of longitudes and latitudes in decimal degrees
        """
        p = self.getproj(withgrid = True)
        lon, lat = p(i + 0.5, j + 0.5, inverse = True)
        return lon, lat
    
    def _newlike(self):
        """
        Internal function to return a file of the same class if a PsueoNetCDFFile
        """
        if isinstance(self, PseudoNetCDFFile):
            outt = type(self)
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        return outf
        
    def renameVariable(self, oldkey, newkey, inplace = False, copyall = True):
        """
        Rename variable (oldkey)
        
        Parameters
        ----------
        oldkey : variable to be renamed
        newkey : new dame for variable
        inplace : create the new variable in this netcdf file (default False)
        copyall : if not inplace, should all variables be copied to new file
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with renamed variable (this file if inplace = True)
        """
        return self.renameVariables(**{oldkey: newkey})
    
    def renameVariables(self, inplace = False, copyall = True, **newkeys):
        """
        Rename variables for each oldkey: newkey dictionary item
        
        Parameters
        ----------
        newkeys : dictionary where key is the oldkey and varlue is the newkey
        inplace : create the new variable in this netcdf file (default False)
        copyall : if not inplace, should all variables be copied to new file
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with renamed variable (this file if inplace = True)
        """
        if inplace:
            outf = self
        else:
            outf = self._copywith(props = True, dimensions = True, variables = copyall, data = copyall)
        
        for oldkey, newkey in newkeys.items():
            outf.copyVariable(self.variables[oldkey], key = newkey)
            if oldkey in outf.variables:
                del outf.variables[oldkey]
        
        return outf
            
    def renameDimension(self, oldkey, newkey, inplace = False):
        """
        Rename dimension (oldkey) in dimensions and in all variables
        
        Parameters
        ----------
        oldkey : dimension to be renamed
        newkey : new dame for dimension
        inplace : create the new variable in this netcdf file (default False)
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with renamed variable (this file if inplace = True)
        """
        return self.renameDimensions(**{oldkey: newkey})
    
    def renameDimensions(self, inplace = False, **newkeys):
        """
        Rename dimension (oldkey) in dimensions and in all variables
        
        Parameters
        ----------
        newkeys : dictionary where key is the oldkey and varlue is the newkey
        inplace : create the new variable in this netcdf file (default False)
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with renamed variable (this file if inplace = True)
        """
        if inplace:
            outf = self
        else:
            outf = self.copy()
        
        for oldkey, newkey in newkeys.items():
            outf.dimensions[newkey] = outf.dimensions[oldkey]

        for k, v in outf.variables.items():
            olddims = v.dimensions
            newdims = tuple([newkeys.get(dk, dk) for dk in olddims])
            if newdims != olddims:
                v.dimensions = newdims
        for oldkey, newkey in newkeys.items():
            del outf.dimensions[oldkey]
        
        return outf
    
    def eval(self, expr, inplace = False, copyall = False):
        """
        Evaluate expr and return a PseudoNetCDFFile object with resutl
        
        Parameters
        ----------
        expr : string with expression to evaluate
        inplace : create the new variable in this netcdf file (default False)
        copyall : if not inplace, should all variables be copied to new file
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with renamed variable (this file if inplace = True)
        """
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
            if copyall: newkeys = None
            else: newkeys = [key]
            outf = self.subsetVariables(newkeys)
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
        """
        Set ncattrs from attdict keys and values
        
        Parameters
        ----------
        attdict : dictionary of properties
        
        Returns
        -------
        None
        """
        for pk, pv in attdict.items():
            setattr(self, pk, pv)
    
    def getncatts(self):
        """
        Return all ncattrs keys and values as a dictionary

        Returns
        ----------
        attdict : dictionary of properties
        """
        outd = OrderedDict()
        for pk in self.ncattrs():
            outd[pk] = getattr(self, pk)
        return outd

    def _copywith(self, props = True, dimensions = True, variables = False, data = False):
        """
        Internal function for making copies of the same type
        
        Parameters
        ----------
        props : boolean include properties (default: True)
        dimensions : boolean include dimensions (default: True)
        variables : boolean include variable structures (default: False)
        data : boolean include variable data (default: False)
        
        Returns
        -------
        outf : PseudoNetCDFFile instance
        
        Notes
        -----
        Internal function does not return variables by default.
        This is useful for functions like slice, apply, eval, etc.
        
        The _ in _copywith means this is a private function and the 
        call interface may change.
        """
        outf = self._newlike()
        if props:
            for pk in self.ncattrs():
                setattr(outf, pk, getattr(self, pk))
        if dimensions:
            for dk, dv in self.dimensions.items():
                newdl = len(dv)
                ndv = outf.createDimension(dk, newdl) 
                ndv.setunlimited(dv.isunlimited())
        if variables: 
            for vk, vv in self.variables.items():
                outf.copyVariable(vv, key = vk, withdata = data)
        return outf
    
    def copy(self, props = True, dimensions = True, variables = True, data = True):
        """
        Function for making copies of the same type
        
        Parameters
        ----------
        props : boolean include properties (default: True)
        dimensions : boolean include dimensions (default: True)
        variables : boolean include variable structures (default: True)
        data : boolean include variable data (default: True)
        
        Returns
        -------
        outf : PseudoNetCDFFile instance
        """
        return self._copywith(props = props, dimensions = dimensions, variables = variables, data = data)
        
    def applyAlongDimensions(self, **dimfuncs):
        """
        Similar to numpy.apply_along_axis, but for damed dimensions and 
        processes dimensions as well as variables
        
        Parameters
        ----------
        dimfuncs : key value pairs where the key is a dimenions and the value
                   is a 1D function (func1d) or a dictionary. If the value is a dictionary
                   it must include func1d as a function and any keyword arguments
                   as additional options
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with variables and dimensions after processing
        """
        outf = self._copywith(props = True, dimensions = True)
        for dk, df in dimfuncs.items():
            dv = self.dimensions[dk]
            if dk in dimfuncs:
                if dk in self.variables:
                    dvar = self.variables[dk]
                else:
                    dvar = np.arange(len(dv))
                if isinstance(df, str):
                    newdl = getattr(dvar[...], df)(keepdims = True).size
                else: newdl = df(dvar[:]).size
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
                         noopts = False
                     else:
                         opts['func1d'] = dfunc
                         noopts = True
                     if noopts and isinstance(dfunc, str):
                         newvals = getattr(newvals, dfunc)(axis = di, keepdims = True)
                     else:
                         newvals = np.apply_along_axis(dfunc, di, newvals)
             newvaro = outf.copyVariable(varo, key = vark, withdata = False)
             newvaro[...] = newvals
        
        return outf
    
    def getTimes(self):
        """
        Get an array of datetime objects
        
        Notes
        -----
        self must have a time or TFLAG variable
        """
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
    
    def stack(self, other, stackdim):
        """
        Concatenates all variables on stackdim
        
        Parameters
        ----------
        other : netcdf-like object
        stackdim : dimension name
        
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with stacked variables and dimension equal to new lenght
        """
        outf = self._copywith(props = True, dimensions = False)
        fs = [self, other]
        dimensions = [f_.dimensions for f_ in fs]
        shareddims = {}
        for dimk, dim in self.dimensions.items():
            if dimk == stackdim:
                continue
            dimlens = [len(dims[dimk]) for dims in dimensions]
            if all([len(dim) == i for i in dimlens]):
                shareddims[dimk] = len(dim)
        differentdims = [set(dims.keys()).difference(shareddims.keys()) for dims in dimensions]
        assert(all([different.union([stackdim]) == set([stackdim]) for different in differentdims]))
        for dimkey in shareddims:
            ind = self.dimensions[dimkey]
            outd = outf.createDimension(dimkey, len(ind))
            outd.setunlimited(ind.isunlimited())
        outd = outf.createDimension(stackdim, sum([len(dims[stackdim]) for dims in dimensions]))
        outd.setunlimited(self.dimensions[stackdim].isunlimited())
        for tmpf in fs:
            for varkey, var in tmpf.variables.items():
                if not stackdim in var.dimensions:
                    if varkey in outf.variables:
                        if np.array_equal(outf.variables[varkey][...], var[...]):
                            pass
                        elif not varkey in self.dimensions:
                            warn('Got duplicate variables for %s without stackable dimension; first value retained' % varkey)
                        continue
                    else:
                        outvals = var[...]
                else:
                    if not varkey in outf.variables.keys():
                        axisi = list(var.dimensions).index(stackdim)
                        outvals = np.ma.concatenate([f_.variables[varkey][:] for f_ in fs], axis = axisi)
                    else: continue
                outvar = outf.copyVariable(var, key = varkey, withdata = False)
                outvar[...] = outvals
        
        return outf
     
    def subsetVariables(self, varkeys, inplace = False, exclude = False):
        """
        Return a PseudoNetCDFFile with only varkeys
        
        Parameters
        ----------
        varkeys : iterable of keys to keep
        inplace : if true (default false), then remove other variable from this file
        exclude : if True (default False), then remove just these variables
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with variables
        """
        if exclude:
            varkeys = list(set(list(self.variables)).difference(varkeys))
        if inplace:
            outf = self
            for varkey in list(outf.variables):
                if not varkey in varkeys:
                    del outf.variables[varkey]
        else:
            outf = self._copywith(props = True, dimensions = True)
            for varkey in varkeys:
                varo = self.variables[varkey]
                newvaro = outf.copyVariable(varo, key = varkey, withdata = False)
                newvaro[...] = varo[...]
        return outf 

    def sliceDimensions(self, newdims = ('POINTS',), **dimslices):
        """
        Return a netcdflike object with dimensions sliced
        
        Parameters
        ----------
        dimslices : key value pairs where the key is a dimension and the
                    value is a valid slice object (slices, ints or iterables)
                    if iterables are provided, all iterables must be the same
                    size and shape. If the arrays are not 1D, newdims must have ndim
                    names
        newdims : names for new dimensions. When more than one iterable applies to
                 a variable slice, fancy indexing removes both dimensions and creates
                 a new one of the iterable lengths
        Returns
        -------
        outf : PseudoNetCDFFile instance with variables and dimensions sliced
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
             if anyisarray and needsfancy:
                 concatax = np.argmax(isdarray)
                 odims = [dk for dk in vdims if not isarray.get(dk, False)]
                 for newdim in newdims[::-1]:
                     odims.insert(concatax, newdim)
             
             newvaro = outf.copyVariable(varo, key = vark, dimensions = odims, withdata = False)
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
             try: newvaro[...] = newvals
             except: newvaro[...] = newvals.reshape(newvaro.shape)
        
        return outf
    
    def removeSingleton(self, dimkey = None):
        """
        Return a netcdflike object with dimensions sliced
        
        Parameters
        ----------
        dimkey : key of dimension to be evaluated for removal; if None, evaluate all.
                 only singleton dimensions will be removed.
        
        Returns
        -------
        outf : PseudoNetCDFFile instance with dimensions removed
        """
        outf = self._copywith(props = True, dimensions = False)
        removed_dims = []
        for dk, d in self.dimensions.items():
            ni = len(d)
            if (dimkey is None or dk == dimkey) and ni == 1:
                removed_dims.append(dk)
            else:
                tempd = outf.createDimension(dk, ni)
                tempd.setunlimited(d.isunlimited())

        for vk, v in self.variables.items():
            olddims = v.dimensions
            newdims = tuple([dk for dk in v.dimensions if not dk in removed_dims])
            sdims = tuple([(di, dk) for di, dk in enumerate(olddims) if dk not in newdims])[::-1]
            propd = dict([(pk, getattr(v, pk)) for pk in v.ncattrs()])
            ov = outf.createVariable(vk, v.dtype.char, newdims, **propd)
            outvals = v[...]
            for di, dk in sdims:
                outvals = outvals.take(0, axis = di)

            ov[...] = outvals[...]
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
        """
        True if this file or object can be identified
        for use by this class. Useful to override for
        classes that can be initialized from disk.
        """
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
        Create a dimension
        
        Parameters
        ----------
        name : string name for dimension
        length : maximum length of dimension

        Returns
        -------
        dim : new dimension
        """
        dim = self.dimensions[name] = PseudoNetCDFDimension(self, name, length)
        return dim
    
    def copyVariable(self, var, key = None, dtype = None, dimensions = None, fill_value = None, withdata = True):
        """
        Copy var into self as vark
        
        Parameters
        ----------
        var : netCDF4.Variable-like object (must have ncattrs and setncatts)
        key : key for variable in self (can be omitted if var has name,
              standard_name, or long_name)
        dtype : change the data type to dtype
        dimensions : change the dimensions to dimensions
        fill_value : change the fill_value to
        withdata : default True, copies data
        
        Returns
        -------
        myvar : copy of var
        """
        if key is None:
            for propk in ['name', 'standard_name', 'long_name']:
                if hasattr(var, propk):
                    key = getattr(var, propk)
            else:
                raise AttributeError('varkey must be supplied because var has no name, standard_name or long_name')
        
        if dtype is None:
            dtype = var.dtype
        if dimensions is None:
            dimensions = var.dimensions
        if fill_value is None:
            for pk in ('fill_value', 'missing_value', '_FillValue'):
                fill_value = getattr(var, pk, None)
                if not fill_value is None: break
        
        myvar = self.createVariable(key, dtype, dimensions, fill_value = fill_value)
        attrs = OrderedDict()
        for propk in var.ncattrs():
            attrs[propk] = getattr(var, propk)
        myvar.setncatts(attrs)
        if withdata: myvar[:] = var[:]
        return myvar
     
    def createVariable(self, name, type, dimensions, fill_value = None, **properties):
        """
        Create a variable
        
        Parameters
        ----------
        name : string
        type : numpy dtype code (e.g., 'f', 'i', 'd')
        dimensions : tuple of dimension keys that can be
                     found in objects' dimensions dictionary
        
        Returns
        -------
        var : new variable
        """
        import numpy as np
        
        if fill_value is None:
            for pk in 'missing_value _FillValue'.split():
                fill_value = properties.get(pk, None)
                if not fill_value is None: break
        if type == 'S': type = 'c'
        if isinstance(properties.get('values', 1), np.ma.MaskedArray) or not fill_value is None:
            var = self.variables[name] = PseudoNetCDFMaskedVariable(self, name, type, dimensions, **properties)
        else:
            var = self.variables[name] = PseudoNetCDFVariable(self, name, type, dimensions, **properties)
        return var

    def close(self):
        """
        Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
        """
        pass
    
    def save(self, *args, **kwds):
        """
        Provides access to pncwrite for self
        
        Parameters
        ----------
        see Help pncwrite

        Returns
        -------
        see Help pncwrite
        """
        from PseudoNetCDF import pncwrite
        return pncwrite(self, *args, **kwds)
    
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

class netcdf(PseudoNetCDFFile, NetCDFFile):
    def createDimension(self, *args, **kwds):
        return NetCDFFile.createDimension(self, *args, **kwds)
    
    def createVariable(self, *args, **kwds):
        return NetCDFFile.createVariable(self, *args, **kwds)
    
    def __setattr__(self, k, v):
        NetCDFFile.__setattr__(self, k, v)
        
    def __delattr__(self, k):
        NetCDFFile.__delattr__(self, k)
    
    def __new__(cls, *args, **kwds):
        return NetCDFFile.__new__(cls, *args, **kwds)
    
    def __init__(self, *args, **kwds):
        NetCDFFile.__init__(self, *args, **kwds)
    
    def ncattrs(self):
        return NetCDFFile.ncattrs(self)
    
    def _newlike(self):
        """
        Internal function to return a file of the same class if a PsueoNetCDFFile
        """
        outf = PseudoNetCDFFile()
        return outf
    
    @classmethod
    def isMine(cls, path, *args, **kwds):
        """
        True if this file or object can be identified
        for use by this class. Useful to override for
        classes that can be initialized from disk.
        """
        tmpf = open(path, 'rb')
        tmpf.seek(0,0)
        cdftest = tmpf.read(3)
        tmpf.seek(1, 0)
        hdftest = tmpf.read(3)
        tmpf.close()
        if cdftest == b'CDF':
            return True
        elif hdftest == b'HDF':
            try:
                f = cls(path, *args, **kwds)
                return True
            except:
                return False
        else:
            return False


registerreader('nc', netcdf)
registerreader('ncf', netcdf)

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
    
    def _makencf(self):
        from numpy import arange
        tncf = self.tncf  = PseudoNetCDFFile()

        tncf.createDimension('TIME', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        o3 = tncf.createVariable('O3', 'f', ('TIME', 'LAY', 'ROW', 'COL'))

        o3[:] = arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)
        o3.units = 'ppbv'
        return tncf

    def testCopyVariable(self):
        tncf = self._makencf()
        var = tncf.copyVariable(tncf.variables['O3'], key = 'O3_PPB', withdata = True)
        self.assertEqual(True, (var[:] == tncf.variables['O3_PPB']).all())
        
    def testSubsetVariables(self):
        tncf = self._makencf()
        var = tncf.copyVariable(tncf.variables['O3'], key = 'O3_PPB', withdata = True)
        var[:] *= 1e3
        var = tncf.copyVariable(tncf.variables['O3'], key = 'O3_PPT', withdata = True)
        var[:] *= 1e6
        sncf = tncf.subsetVariables(['O3_PPT'])
        self.assertEqual(len(sncf.variables), 1)
        self.assertEqual(set(sncf.variables), set(['O3_PPT']))
        sncf = tncf.subsetVariables(['O3_PPT'], exclude = True)
        self.assertEqual(len(sncf.variables), 2)
        self.assertEqual(set(sncf.variables), set(['O3', 'O3_PPB']))
        # I'm here BHH

    def testRenameVariables(self):
        tncf = self._makencf()
        sncf = tncf.renameVariables(O3 = 'O3_PPM')
        self.assertEqual(len(sncf.variables), 1)
        self.assertEqual(set(sncf.variables), set(['O3_PPM']))
        
    def testRenameDimensions(self):
        tncf = self._makencf()
        sncf = tncf.renameDimensions(TIME = 'TSTEP')
        self.assertEqual(len(sncf.dimensions), len(tncf.dimensions))
        self.assertEqual(set(sncf.dimensions), set(['TSTEP', 'LAY', 'ROW', 'COL']))
        
    def testSliceDimension(self):
        tncf = self._makencf()
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(TIME = 0)
        self.assertEqual(len(sncf.dimensions['TIME']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == tncf.variables['O3'][0]).all())
        sncf = tncf.sliceDimensions(TIME = [0])
        self.assertEqual(len(sncf.dimensions['TIME']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == tncf.variables['O3'][0]).all())
        sncf = tncf.sliceDimensions(TIME = [0, 8])
        self.assertEqual(len(sncf.dimensions['TIME']), 2)
        self.assertEqual(True, (sncf.variables['O3'][:] == tncf.variables['O3'][[0, 8]]).all())
        sncf = tncf.sliceDimensions(TIME = [0, 8], ROW = 2, COL = 3)
        self.assertEqual(len(sncf.dimensions['TIME']), 2)
        self.assertEqual(len(sncf.dimensions['ROW']), 1)
        self.assertEqual(len(sncf.dimensions['COL']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == o3[[0, 8], :, 2, 3][:, :, None, None]).all())
        i = np.arange(4)
        sncf = tncf.sliceDimensions(TIME = i, LAY = i, ROW = i, COL = i)
        self.assertEqual(len(sncf.dimensions['POINTS']), 4)
        self.assertEqual(True, (sncf.variables['O3'][:] == o3[i, i, i, i]).all())
        
    def testApplyAlongDimensions(self):
        tncf = self._makencf()
        o3 = tncf.variables['O3'][:]
        ancf = tncf.applyAlongDimensions(LAY = 'min')
        self.assertEqual(True, (ancf.variables['O3'][:] == o3.min(1, keepdims = True)).all())
        # Testing convolution; useful for mda8
        ancf = tncf.applyAlongDimensions(TIME = lambda x: np.convolve(x, np.ones(2, dtype = 'f') / 2., mode = 'valid'))
        co3 = (o3[1:] + o3[:-1]) / 2
        self.assertEqual(True, (ancf.variables['O3'][:] == co3).all())
        ancf = tncf.applyAlongDimensions(TIME = lambda x: np.convolve(x, np.ones(2, dtype = 'f') / 2., mode = 'valid')).applyAlongDimensions(TIME = np.max)
        mco3 = co3.max(0, keepdims = True)
        self.assertEqual(True, (ancf.variables['O3'][:] == mco3).all())
        
        
    def testNetCDFFileNew(self):
        t = PseudoNetCDFFile.__new__(PseudoNetCDFFile)
        self.assertEqual(t.variables, {})
        self.assertEqual(t.dimensions, {})
        self.assertEqual(t.ncattrs(), ())
        
    def testNetCDFFileInit(self):
        from numpy import arange
        self._makencf()
        tncf = self.tncf
        self.assertEqual(len(tncf.dimensions['TIME']), 24)
        self.assertEqual(len(tncf.dimensions['LAY']), 4)
        self.assertEqual(len(tncf.dimensions['ROW']), 5)
        self.assertEqual(len(tncf.dimensions['COL']), 6)

        
        tncf.fish = 2
        setattr(tncf, 'FROG-DOG', 'HAPPY')

        self.assertEqual(set(tncf.variables.keys()), set(['O3']))
        o3 = tncf.variables['O3']
        self.assertEqual(True, (o3 == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        self.assertEqual(o3.typecode(), 'f')

        filedims = list(tncf.dimensions)
        filedims.sort()
        vardims = list(o3.dimensions)
        vardims.sort()

        self.assertEqual(filedims, vardims)
        from PseudoNetCDF.pncgen import Pseudo2NetCDF
        n = Pseudo2NetCDF().convert(tncf)
        self.assertEqual(set(n.variables.keys()), set(['O3']))
        self.assertEqual(dict([(k, len(v)) for k, v in n.dimensions.items()]), {'TIME': 24, 'LAY': 4, 'ROW': 5, 'COL': 6})
        self.assertEqual(True, (n.variables['O3'][...] == tncf.variables['O3'][...]).all())
        self.assertEqual(n.variables['O3'].units, 'ppbv')
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
        self.assertEqual(True, (tncf.variables['O3'] == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        
    def runTest(self):
        pass
    

if __name__ == '__main__':
    unittest.main()
