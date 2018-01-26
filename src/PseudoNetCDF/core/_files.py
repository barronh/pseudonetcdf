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
                lonb = self.variables['longitude_bounds']
                latb = self.variables['latitude_bounds']

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
        newkeys : dictionary where key is the oldkey and value is the newkey
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
        newkeys : dictionary where key is the oldkey and value is the newkey
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
    
    def insertDimension(self, newonly = True, multionly = False, before = None, inplace = False, **newdims):
        """
        Insert dimensions with keys and lengths from newdims
        
        
        Parameters
        ----------
        newdims   : dictionary where key is the new dimension and value is the length
        newonly   : Only add dimension to variables that do not already have it,
                    default True
        multionly : Only add dimension if there are already more than one (good for 
                    ignoring coordinate dimensions)
        before    : if variable has this dimension, insert the new dimension before 
                    it. Otherwise, add to the beginning.
                        
        inplace   : create the new variable in this netcdf file (default False)
        
        Returns
        -------
        outf : PseudoNetCDFFile instance will new dimension in dimensions and variables
        
        Notes
        -----
        
        1. Adding a non unity dimension will cause the data to be repeated
           along the new axis.
        2. If order of addition matters, use multiple calls. newdimsuse
           will be a non-ordered dictionary
        """
        if inplace:
            outf = self
        else:
            outf = self.copy(variables = False)
        
        for dk, dv in newdims.items():
            if not dk in outf.dimensions:
                ndv = outf.createDimension(dk, dv)
                
            for vk, vv in self.variables.items():
                vdims = list(vv.dimensions)
                if (newonly and dk in vdims) or (multionly and len(vdims) == 1):
                    outf.copyVariable(vv, key = vk, withdata = True)
                    continue
                if before is None:
                    ndims = tuple([dk] + vdims)
                    bi = 0
                else:
                    ndims = vdims
                    if before in vdims:
                        bi = vdims.index(before)
                    else:
                        outf.copyVariable(vv, key = vk, withdata = True)
                        continue
                    ndims.insert(bi, dk)
                var = outf.copyVariable(vv, key = vk, dimensions = ndims, withdata = False)
                
                var[...] = np.expand_dims(self.variables[vk][...], axis = bi)
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
        dimfuncs : key value pairs where the key is a dimensions and the value
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
            calendar = getattr(time, 'calendar', 'gregorian')
            if 'since' in time.units:
                unit, base = time.units.strip().split(' since ')
                sdate = _parse_ref_date(base)
                daysinyear = {'noleap': 365, '365_day': 365, 'all_leap': 366, '366_day': 366}
                if calendar in daysinyear:
                    yeardays = daysinyear[calendar]
                    warn('Dates shown in standard calendar with leap year, but calculated with %d days' % yeardays)
                    year = sdate.year
                    month = sdate.month
                    day = sdate.day
                    yearincrs = time[:] / {'days': yeardays, 'hours': yeardays*24, 'minutes': yeardays*24*60, 'seconds': yeardays*24*60}[unit]
                    out = np.array([datetime(year + int(yearinc), month, day) + timedelta(days = (yearinc % 1) * yeardays) for yearinc in yearincrs])
                else:
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
        if withdata:
            try: myvar[:] = var[:]
            except: myvar[...] = var[...]
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
        from datetime import datetime, timedelta
        self.testncf = self._makencf()
        self.mymeta = set(['time', 'latitude', 'longitude', 'latitude_bounds', 'longitude_bounds', 'lambert_conformal_conic'])
        self.myvars = self.mymeta.union(['O3'])
        rtime = datetime.strptime('1970-01-01 00:00:00+0000', '%Y-%m-%d %H:%M:%S%z')
        self.mytimes = np.array([rtime + timedelta(hours = i) for i in range(24)])
    def _makencf(self):
        from numpy import arange
        tncf = PseudoNetCDFFile()

        tncf.createDimension('TSTEP', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)
        tncf.createDimension('nv', 4)
        tncf.str_one = '1'
        tncf.int_two = 2
        tncf.float_threeptfive = 3.5
        tncf.Conventions = 'CF-1.6'
        o3 = tncf.createVariable('O3', 'f', ('TSTEP', 'LAY', 'ROW', 'COL'))

        o3[:] = arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)
        o3.units = 'ppbv'
        o3.grid_mapping = 'lambert_conformal_conic'
        time = tncf.createVariable('time', 'd', ('TSTEP',))
        time.long_name = 'time'
        time.units = 'hours since 1970-01-01 00:00:00+0000'
        time[:] = np.arange(24)
        
        crs = tncf.createVariable('lambert_conformal_conic', 'i', ())
        crs.grid_mapping_name = 'lambert_conformal_conic'
        crs.standard_parallel = np.array([30., 45.])
        crs.longitude_of_central_meridian = -97.
        crs.latitude_of_projection_origin = 40.
        crs.false_northing = 1620000.
        crs.false_easting = 2412000.
        crs.semi_major_axis = 6371000.
        crs.semi_minor_axis = 6371000.
        lat = tncf.createVariable('latitude', 'f', ('ROW', 'COL'))
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lon = tncf.createVariable('longitude', 'f', ('ROW', 'COL'))
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'
        latb = tncf.createVariable('latitude_bounds', 'f', ('ROW', 'COL', 'nv'))
        latb.long_name = 'latitude_bounds'
        latb.units = 'degrees_north'
        lonb = tncf.createVariable('longitude_bounds', 'f', ('ROW', 'COL', 'nv'))
        lonb.long_name = 'longitude_bounds'
        lonb.units = 'degrees_east'
        lon[:] = [[-120.21161038333193, -120.21160114763147, -120.21159191193058, -120.21158267622918, -120.21157344052737, -120.21156420482505], [-120.21161271536134, -120.21160347966001, -120.21159424395826, -120.21158500825604, -120.21157577255335, -120.21156653685021], [-120.21161504739118, -120.21160581168901, -120.21159657598642, -120.21158734028334, -120.2115781045798, -120.2115688688758], [-120.21161737942151, -120.21160814371851, -120.21159890801503, -120.21158967231109, -120.21158043660672, -120.21157120090189], [-120.21161971145229, -120.21161047574842, -120.21160124004409, -120.21159200433934, -120.21158276863409, -120.21157353292838]] 

        lat[:] = [[22.748507533242535, 22.748509683865187, 22.74851183448702, 22.74851398510802, 22.748516135728206, 22.748518286347593], [22.74851605050742, 22.748518201130356, 22.748520351752475, 22.74852250237377, 22.74852465299425, 22.748526803613903], [22.74852456777256, 22.748526718395773, 22.748528869018187, 22.748531019639763, 22.748533170260536, 22.74853532088048], [22.748533085037966, 22.74853523566145, 22.74853738628417, 22.748539536906023, 22.7485416875271, 22.748543838147327], [22.74854160230359, 22.74854375292739, 22.748545903550376, 22.748548054172538, 22.748550204793883, 22.7485523554144]]
        lonb[:] = [[[-120.21161038333193, -120.21160114763147, -120.21160347966001, -120.21161271536134], [-120.21160114763147, -120.21159191193058, -120.21159424395826, -120.21160347966001], [-120.21159191193058, -120.21158267622918, -120.21158500825604, -120.21159424395826], [-120.21158267622918, -120.21157344052737, -120.21157577255335, -120.21158500825604], [-120.21157344052737, -120.21156420482505, -120.21156653685021, -120.21157577255335], [-120.21156420482505, -120.2115549691223, -120.2115573011466, -120.21156653685021]], [[-120.21161271536134, -120.21160347966001, -120.21160581168901, -120.21161504739118], [-120.21160347966001, -120.21159424395826, -120.21159657598642, -120.21160581168901], [-120.21159424395826, -120.21158500825604, -120.21158734028334, -120.21159657598642], [-120.21158500825604, -120.21157577255335, -120.2115781045798, -120.21158734028334], [-120.21157577255335, -120.21156653685021, -120.2115688688758, -120.2115781045798], [-120.21156653685021, -120.2115573011466, -120.21155963317135, -120.2115688688758]], [[-120.21161504739118, -120.21160581168901, -120.21160814371851, -120.21161737942151], [-120.21160581168901, -120.21159657598642, -120.21159890801503, -120.21160814371851], [-120.21159657598642, -120.21158734028334, -120.21158967231109, -120.21159890801503], [-120.21158734028334, -120.2115781045798, -120.21158043660672, -120.21158967231109], [-120.2115781045798, -120.2115688688758, -120.21157120090189, -120.21158043660672], [-120.2115688688758, -120.21155963317135, -120.21156196519657, -120.21157120090189]], [[-120.21161737942151, -120.21160814371851, -120.21161047574842, -120.21161971145229], [-120.21160814371851, -120.21159890801503, -120.21160124004409, -120.21161047574842], [-120.21159890801503, -120.21158967231109, -120.21159200433934, -120.21160124004409], [-120.21158967231109, -120.21158043660672, -120.21158276863409, -120.21159200433934], [-120.21158043660672, -120.21157120090189, -120.21157353292838, -120.21158276863409], [-120.21157120090189, -120.21156196519657, -120.21156429722222, -120.21157353292838]], [[-120.21161971145229, -120.21161047574842, -120.21161280777879, -120.2116220434835], [-120.21161047574842, -120.21160124004409, -120.21160357207363, -120.21161280777879], [-120.21160124004409, -120.21159200433934, -120.21159433636801, -120.21160357207363], [-120.21159200433934, -120.21158276863409, -120.21158510066192, -120.21159433636801], [-120.21158276863409, -120.21157353292838, -120.21157586495535, -120.21158510066192], [-120.21157353292838, -120.21156429722222, -120.21156662924835, -120.21157586495535]]]
        latb[:] = [[[22.748507533242535, 22.748509683865187, 22.748518201130356, 22.74851605050742], [22.748509683865187, 22.74851183448702, 22.748520351752475, 22.748518201130356], [22.74851183448702, 22.74851398510802, 22.74852250237377, 22.748520351752475], [22.74851398510802, 22.748516135728206, 22.74852465299425, 22.74852250237377], [22.748516135728206, 22.748518286347593, 22.748526803613903, 22.74852465299425], [22.748518286347593, 22.748520436966125, 22.748528954232754, 22.748526803613903]], [[22.74851605050742, 22.748518201130356, 22.748526718395773, 22.74852456777256], [22.748518201130356, 22.748520351752475, 22.748528869018187, 22.748526718395773], [22.748520351752475, 22.74852250237377, 22.748531019639763, 22.748528869018187], [22.74852250237377, 22.74852465299425, 22.748533170260536, 22.748531019639763], [22.74852465299425, 22.748526803613903, 22.74853532088048, 22.748533170260536], [22.748526803613903, 22.748528954232754, 22.748537471499613, 22.74853532088048]], [[22.74852456777256, 22.748526718395773, 22.74853523566145, 22.748533085037966], [22.748526718395773, 22.748528869018187, 22.74853738628417, 22.74853523566145], [22.748528869018187, 22.748531019639763, 22.748539536906023, 22.74853738628417], [22.748531019639763, 22.748533170260536, 22.7485416875271, 22.748539536906023], [22.748533170260536, 22.74853532088048, 22.748543838147327, 22.7485416875271], [22.74853532088048, 22.748537471499613, 22.748545988766764, 22.748543838147327]], [[22.748533085037966, 22.74853523566145, 22.74854375292739, 22.74854160230359], [22.74853523566145, 22.74853738628417, 22.748545903550376, 22.74854375292739], [22.74853738628417, 22.748539536906023, 22.748548054172538, 22.748545903550376], [22.748539536906023, 22.7485416875271, 22.748550204793883, 22.748548054172538], [22.7485416875271, 22.748543838147327, 22.7485523554144, 22.748550204793883], [22.748543838147327, 22.748545988766764, 22.748554506034104, 22.7485523554144]], [[22.74854160230359, 22.74854375292739, 22.74855227019359, 22.748550119569494], [22.74854375292739, 22.748545903550376, 22.748554420816852, 22.74855227019359], [22.748545903550376, 22.748548054172538, 22.74855657143929, 22.748554420816852], [22.748548054172538, 22.748550204793883, 22.74855872206093, 22.74855657143929], [22.748550204793883, 22.7485523554144, 22.748560872681754, 22.74855872206093], [22.7485523554144, 22.748554506034104, 22.748563023301763, 22.748560872681754]]]
        return tncf

    def testCopyVariable(self):
        tncf = self.testncf
        var = tncf.copyVariable(tncf.variables['O3'], key = 'O3_PPB', withdata = True)
        self.assertEqual(True, (var[:] == tncf.variables['O3_PPB']).all())
        
    def testSubsetVariables(self):
        tncf = self.testncf.copy()
        var = tncf.copyVariable(tncf.variables['O3'], key = 'O3_PPB', withdata = True)
        var[:] *= 1e3
        var = tncf.copyVariable(tncf.variables['O3'], key = 'O3_PPT', withdata = True)
        var[:] *= 1e6
        sncf = tncf.subsetVariables(['O3_PPT'])
        self.assertEqual(len(sncf.variables), 1)
        self.assertEqual(set(sncf.variables), set(['O3_PPT']))
        sncf = tncf.subsetVariables(['O3_PPT'], exclude = True)
        self.assertEqual(len(sncf.variables), len(self.myvars.union(['O3_PPB'])))
        self.assertEqual(set(sncf.variables), self.myvars.union(['O3_PPB']))
        

    def testRenameVariables(self):
        tncf = self.testncf
        sncf = tncf.renameVariables(O3 = 'O3_PPM')
        self.assertEqual(len(sncf.variables), len(self.myvars))
        self.assertEqual(set(sncf.variables), self.mymeta.union(['O3_PPM']))
        
    def testRenameDimensions(self):
        tncf = self.testncf
        sncf = tncf.renameDimensions(TSTEP = 'TIME')
        self.assertEqual(len(sncf.dimensions), len(tncf.dimensions))
        self.assertEqual(set(sncf.dimensions), set(['TIME', 'LAY', 'ROW', 'COL', 'nv']))
        
    def testSliceDimension(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(TSTEP = 0)
        self.assertEqual(len(sncf.dimensions['TSTEP']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == tncf.variables['O3'][0]).all())
        sncf = tncf.sliceDimensions(TSTEP = [0])
        self.assertEqual(len(sncf.dimensions['TSTEP']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == tncf.variables['O3'][0]).all())
        sncf = tncf.sliceDimensions(TSTEP = [0, 8])
        self.assertEqual(len(sncf.dimensions['TSTEP']), 2)
        self.assertEqual(True, (sncf.variables['O3'][:] == tncf.variables['O3'][[0, 8]]).all())
        sncf = tncf.sliceDimensions(TSTEP = [0, 8], ROW = 2, COL = 3)
        self.assertEqual(len(sncf.dimensions['TSTEP']), 2)
        self.assertEqual(len(sncf.dimensions['ROW']), 1)
        self.assertEqual(len(sncf.dimensions['COL']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == o3[[0, 8], :, 2, 3][:, :, None, None]).all())
        i = np.arange(4)
        sncf = tncf.sliceDimensions(TSTEP = i, LAY = i, ROW = i, COL = i)
        self.assertEqual(len(sncf.dimensions['POINTS']), 4)
        self.assertEqual(True, (sncf.variables['O3'][:] == o3[i, i, i, i]).all())
        
    def testApplyAlongDimensions(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        ancf = tncf.applyAlongDimensions(LAY = 'min')
        self.assertEqual(True, (ancf.variables['O3'][:] == o3.min(1, keepdims = True)).all())
        # Testing convolution; useful for mda8
        ancf = tncf.applyAlongDimensions(TSTEP = lambda x: np.convolve(x, np.ones(2, dtype = 'f') / 2., mode = 'valid'))
        co3 = (o3[1:] + o3[:-1]) / 2
        self.assertEqual(True, (ancf.variables['O3'][:] == co3).all())
        ancf = tncf.applyAlongDimensions(TSTEP = lambda x: np.convolve(x, np.ones(2, dtype = 'f') / 2., mode = 'valid')).applyAlongDimensions(TSTEP = np.max)
        mco3 = co3.max(0, keepdims = True)
        self.assertEqual(True, (ancf.variables['O3'][:] == mco3).all())
        
    def testGetMap(self):
        tncf = self.testncf
        m = tncf.getMap(maptype = 'basemap_auto')

    def testGetproj(self):
        tncf = self.testncf
        p = tncf.getproj(withgrid = False, projformat = 'pyproj')

    def testLl2xy(self):
        tncf = self.testncf
        crs = tncf.variables['lambert_conformal_conic']
        y0 = crs.false_northing
        x0 = crs.false_easting
        lon0, lat0 = tncf.xy2ll(x0, y0)
        x0t, y0t = tncf.ll2xy(lon0, lat0)
        self.assertEqual((x0, y0), (x0t, y0t))

    def testLl2ij(self):
        tncf = self.testncf
        crs = tncf.variables['lambert_conformal_conic']
        y0 = 0
        x0 = 0
        lon0, lat0 = tncf.xy2ll(x0, y0)
        i0, j0 = tncf.ll2ij(lon0, lat0)
        self.assertEqual(0, i0)
        self.assertEqual(0, j0)

    def testXy2ll(self):
        tncf = self.testncf
        crs = tncf.variables['lambert_conformal_conic']
        y0 = crs.false_northing
        x0 = crs.false_easting
        lonmid, latmid = tncf.xy2ll(x0, y0)
        self.assertEqual(True, np.allclose(lonmid, crs.longitude_of_central_meridian))
        self.assertEqual(True, np.allclose(latmid, crs.latitude_of_projection_origin))
        lon0, lat0 = tncf.xy2ll(0, 0)
        self.assertEqual(True, np.allclose(lon0, tncf.variables['longitude_bounds'][0, 0, 0]))
        self.assertEqual(True, np.allclose(lat0, tncf.variables['latitude_bounds'][0, 0, 0]))
        
    def testIj2ll(self):
        tncf = self.testncf
        lon0, lat0 = tncf.ij2ll(0, 0)
        self.assertEqual(True, np.allclose(lon0, tncf.variables['longitude'][0, 0]))
        self.assertEqual(True, np.allclose(lat0, tncf.variables['latitude'][0, 0]))

    def testEval(self):
        tncf = self.testncf.copy()
        tncf.eval('O3_PPB = O3 * 1000.', inplace = True, copyall = False)
        o3ppmv = tncf.variables['O3']
        o3ppbv = tncf.variables['O3_PPB']
        self.assertEqual(True, ((o3ppmv * 1000.) == o3ppbv).all())
        

    def testSetncatts(self):
        tncf = self.testncf.copy()
        tncf.setncatts({'test_new1': 1, 'test_new2': 'five'})
        self.assertEqual(tncf.test_new1, 1)
        self.assertEqual(tncf.test_new2, 'five')

    def testGetncatts(self):
        tncf = self.testncf
        t = tncf.getncatts()
        self.assertEqual(t['str_one'], '1')
        self.assertEqual(t['int_two'], 2)
        self.assertEqual(t['float_threeptfive'], 3.5)
        

    def testCopy(self):
        tncf = self.testncf
        nncf = tncf.copy(props = True, dimensions = True, variables = True, data = True)
        for pk in tncf.ncattrs():
            self.assertEqual(getattr(tncf, pk), getattr(nncf, pk, None))
        
        for dk, dv in tncf.dimensions.items():
            dvl = len(dv)
            ndvl = len(nncf.dimensions[dk])
            self.assertEqual(dvl, ndvl)
        
        for vk, vv in tncf.variables.items():
            nvv = nncf.variables[vk]
            self.assertEqual(True, np.allclose(vv[...], nvv[...]))
            self.assertEqual(vv.dimensions, nvv.dimensions)
            self.assertEqual(vv.dtype.char, nvv.dtype.char)
            for pk in vv.ncattrs():
                pv = getattr(vv, pk)
                npv = getattr(nvv, pk, None)
                testv = pv == npv
                self.assertEqual(True, np.any(testv))
            

    def testGetTimes(self):
        tncf = self.testncf
        t = tncf.getTimes()
        self.assertEqual(True, (t == self.mytimes).all())

    def testStack(self):
        tncf = self.testncf
        sncf = tncf.stack(tncf, 'TSTEP')
        to3 = tncf.variables['O3'][:]
        no3 = sncf.variables['O3'][:]
        origlen = len(tncf.dimensions['TSTEP'])
        self.assertEqual(origlen*2, len(sncf.dimensions['TSTEP']))
        self.assertEqual(True, np.allclose(to3, no3[:origlen]))
        self.assertEqual(True, np.allclose(to3, no3[origlen:]))

    def testRemoveSingleton(self):
        tncf = self.testncf
        nncf = tncf.copy()
        nncf.createDimension('test', 1)
        nncf = nncf.removeSingleton(dimkey = 'test')
        self.assertEqual(set(tncf.dimensions), set(nncf.dimensions))
        nncf.createDimension('test', 1)
        nncf = nncf.removeSingleton(dimkey = None)
        self.assertEqual(set(tncf.dimensions), set(nncf.dimensions))

    def testCreateDimension(self):
        tncf = self.testncf.copy()
        ndims = len(tncf.dimensions)
        olddims = set(tncf.dimensions)
        newdims = olddims.union(['newd'])
        tncf.createDimension('newd', 5)
        self.assertEqual(len(tncf.dimensions), ndims + 1)
        self.assertEqual(set(tncf.dimensions), newdims)
        
    def testCreateVariable(self):
        tncf = self.testncf.copy()
        nvars = len(tncf.variables)
        var = tncf.createVariable('test', 'f', ('TSTEP',))
        self.assertEqual(len(tncf.variables), nvars + 1)
        self.assertEqual(var.dimensions, ('TSTEP',))
        self.assertEqual(var.dtype.char, 'f')

    def testSave(self):
        pass
        

    def testNcattrs(self):
        tncf = self.testncf
        tncf.ncattrs()

    def testSetncattr(self):
        tncf = self.testncf
        tncf.setncattr('test', 1)
        self.assertEqual(tncf.test, 1)

    def testDelncattr(self):
        tncf = self.testncf.copy()
        tncf.setncattr('test', 1)
        tncf.delncattr('test')
        self.assertEqual(False, hasattr(tncf, 'test'))
    
    def testInsertDimension(self):
        tncf = self.testncf.copy()
        nncf = tncf.insertDimension(TEST = 2, before = 'LAY')
        self.assertEqual(True, (nncf.variables['O3'][:,0] == nncf.variables['O3'][:,1]).all())
        nncf = tncf.insertDimension(TEST = 2)
        self.assertEqual(True, (nncf.variables['O3'][0,:] == nncf.variables['O3'][1,:]).all())
        
        
        
    def testNetCDFFileNew(self):
        t = PseudoNetCDFFile.__new__(PseudoNetCDFFile)
        self.assertEqual(t.variables, {})
        self.assertEqual(t.dimensions, {})
        self.assertEqual(t.ncattrs(), ())
        
    def testNetCDFFileInit(self):
        from numpy import arange
        self._makencf()
        tncf = self.testncf
        self.assertEqual(len(tncf.dimensions['TSTEP']), 24)
        self.assertEqual(len(tncf.dimensions['LAY']), 4)
        self.assertEqual(len(tncf.dimensions['ROW']), 5)
        self.assertEqual(len(tncf.dimensions['COL']), 6)
        self.assertEqual(len(tncf.dimensions['nv']), 4)

        
        tncf.fish = 2
        setattr(tncf, 'FROG-DOG', 'HAPPY')

        self.assertEqual(set(tncf.variables.keys()), self.myvars)
        o3 = tncf.variables['O3']
        self.assertEqual(True, (o3 == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        self.assertEqual(o3.typecode(), 'f')

        filedims = list(tncf.dimensions)
        filedims.sort()
        vardims = list(o3.dimensions)
        vardims.sort()
        filedims.remove('nv')

        self.assertEqual(filedims, vardims)
        from PseudoNetCDF.pncgen import Pseudo2NetCDF
        n = Pseudo2NetCDF().convert(tncf)
        self.assertEqual(set(n.variables.keys()), self.myvars)
        self.assertEqual(dict([(k, len(v)) for k, v in n.dimensions.items()]), {'TSTEP': 24, 'LAY': 4, 'ROW': 5, 'COL': 6, 'nv': 4})
        self.assertEqual(True, (n.variables['O3'][...] == tncf.variables['O3'][...]).all())
        self.assertEqual(n.variables['O3'].units, 'ppbv')
        self.assertEqual(n.fish, 2)
        self.assertEqual(getattr(n, 'FROG-DOG'), 'HAPPY')
        
    def testNetCDFVariables(self):
        from numpy import arange
        tncf = PseudoNetCDFFile()
        tncf.createDimension('TSTEP', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        const = lambda *args, **kwds: PseudoNetCDFVariable(tncf, args[0], 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), values = arange(24 * 4 * 5 * 6).reshape((24, 4, 5, 6)))
        tncf.variables = PseudoNetCDFVariables(const, ['NO', 'O3'])
        self.assertEqual(True, (tncf.variables['O3'] == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        
    def runTest(self):
        pass
    

if __name__ == '__main__':
    unittest.main()
