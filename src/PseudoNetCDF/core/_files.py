import unittest
from PseudoNetCDF._getreader import registerreader
from PseudoNetCDF.netcdf import NetCDFFile, NetCDFVariable
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

class PseudoNetCDFFile(PseudoNetCDFSelfReg, object):
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
                if lonb.ndim == 3:
                    llcrnrlon = lonb[0, 0, 0]
                    urcrnrlon = lonb[-1, -1, 2]
                elif lonb.ndim == 2:
                    llcrnrlon = lonb[0, 0]
                    urcrnrlon = lonb[-1, -1]
                elif lonb.ndim == 1:
                    llcrnrlon = lonb[0]
                    urcrnrlon = lonb[-1]
                if latb.ndim == 3:
                    llcrnrlat = latb[0, 0, 0]
                    urcrnrlat = latb[-1, -1, 2]
                elif latb.ndim == 2:
                    llcrnrlat = latb[0, 0]
                    urcrnrlat = latb[-1, -1]
                elif latb.ndim == 1:
                    llcrnrlat = latb[0]
                    urcrnrlat = latb[-1]
                edges = dict(llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat,
                             urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat)
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
    
    def _getydim(self):
        for dk in 'latitude lat south_north ROW y'.split():
            if dk in self.dimensions:
                 return dk
    
    def _getxdim(self):
        for dk in 'longitude long west_east COL x'.split():
            if dk in self.dimensions:
                 return dk
    
    def ll2ij(self, lon, lat, bounds = 'ignore'):
        """
        Converts lon/lat to 0-based indicies (0,M), (0,N)
        
        Parameters
        ----------
        lon : scalar or iterable of longitudes in decimal degrees
        lat : scalar or iterable of latitudes in decimal degrees
        bounds : ignore, error, warn if i,j are out of domain
        
        Returns
        -------
        i, j : indices (0-based) for variables
        """
        import numpy as np
        p = self.getproj(withgrid = True)
        x, y = p(lon, lat)
        i = np.asarray(x).astype('i')
        j = np.asarray(y).astype('i')
        if bounds == 'ignore':
            pass
        else:
            nx = len(self.dimensions[self._getxdim()])
            ny = len(self.dimensions[self._getydim()])
            lowi = (i < 0)
            lowj = (j < 0)
            highi = (i >= nx)
            highj = (j >= ny)
            outb = (lowi | lowj | highi | highj)
            nout = outb.sum()
            if nout > 0:
                message ='{} Points out of bounds; {}'.format(nout, np.where(outb))
                if bounds == 'error':
                    raise ValueError(message)
                else:
                    pncwarn(message)
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
    
    def insertDimension(self, newonly = True, multionly = False, before = None, after = None, inplace = False, **newdims):
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
                    it. Otherwise, add to the beginning. (before take precedence

        after     : if variable has this dimension, insert the new dimension after 
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
                ndims = [_dk for _dk in vdims]
                if before in vdims:
                    bi = vdims.index(before)
                elif after in vdims:
                    bi = vdims.index(after) + 1
                elif not (before is None and after is None):
                        outf.copyVariable(vv, key = vk, withdata = True)
                        continue
                else:
                    bi = 0
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
                if isinstance(tmpvar, (PseudoNetCDFVariable, NetCDFVariable)):
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
        propd['expression'] = expr
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
    
    def plot(self, varkey, plottype = 'longitude-latitude', ax_kw = {}, plot_kw = {}, cbar_kw = {}, dimreduction = 'mean'):
        """
        Parameters
        ----------
        self : the PseudoNetCDF file instance
        varkey : the variable to plot
        plottype : longitude-latitude, latitude-pressure, longitude-pressure, vertical-profile,
                   time-longitude, time-latitude, time-pressure, default, longitude-latitude
        ax_kw : keywords for the axes to be created
        plot_kw : keywords for the plot (plot, scatter, or pcolormesh) to be created
        cbar_kw : keywords for the colorbar
        """
        import matplotlib.pyplot as plt
        from ..coordutil import getbounds
        apply2dim = {}
        var = self.variables[varkey]
        varunit = varkey
        if hasattr(var, 'units'):
            varunit += var.units.strip()
        
        dimlens = dict([(dk, len(self.dimensions[dk])) for dk in var.dimensions])
        dimpos = dict([(dk, di) for di, dk in enumerate(var.dimensions)])
        xkey, ykey = plottype.split('-')
        if not ykey == 'profile':
            for dimkey in list(dimlens):
                if not dimkey in (xkey, ykey) and dimlens.get(dimkey, 1) > 1:
                    apply2dim[dimkey] = dimreduction
        
        if len(apply2dim) > 0:
           myf = self.applyAlongDimensions(**apply2dim)
           var = myf.variables[varkey]
           dimlens = dict([(dk, len(self.dimensions[dk])) for dk in var.dimensions])
        else:
           myf = self
        if ykey in ('profile',):
           vaxi = var.dimensions.index(xkey)
           vsize = var.shape[vaxi]
           vals = np.rollaxis(var[:], layaxi).reshape(laysize, -1)
        else:
           vals = var[:].squeeze()
        
        if xkey == 'time':
            xm = myf.getTimes()
            dx = np.diff(xm)[-1]
            x = np.append(xm, xm[-1] + dx)
            x = plt.matplotlib.dates.date2num(x)
        else:
            x = getbounds(myf, xkey)
        
        ax = plt.gca(**ax_kw)
        if ykey in ('profile',):
            y = getbounds(myf, xkey)
            x0 = vals[:].min(0)
            xm = vals[:].mean(0)
            x1 = vals[:].max(0)
            ax.fill_betweenx(y = y, x0 = x0, x1 = x1, label = varkey + '(min, max)')
            ax.plot(xm, y, label = varkey, **plot_kw)
            ax.set_ylabel(xkey)
            ax.set_xlabel(varunit)
            return ax
             
        if ykey == 'time':
            ym = myf.getTimes()
            dy = np.diff(ym)[-1]
            y = np.append(ym, ym[-1] + dy)
            y = plt.matplotlib.dates.date2num(y)
        else:
            y = getbounds(myf, ykey)
        
        if dimpos[xkey] < dimpos[ykey]:
            vals = vals.T
        p = ax.pcolormesh(x, y, vals, **plot_kw)
        ax.figure.colorbar(p, label = varunit, **cbar_kw)
        if xkey == 'time':
            ax.xaxis.set_major_formatter(plt.matplotlib.dates.AutoDateFormatter(plt.matplotlib.dates.AutoDateLocator()))
        if ykey == 'time':
            ax.yaxis.set_major_formatter(plt.matplotlib.dates.AutoDateFormatter(plt.matplotlib.dates.AutoDateLocator()))
        if plottype == 'longitude-latitude':
            try:
                bmap = myf.getMap()
                bmap.drawcoastlines(ax = ax)
                bmap.drawcountries(ax = ax)
            except:
                pass
        else:
            ax.set_xlabel(xkey)
            ax.set_ylabel(ykey)
        return ax
    
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
    
    def interpDimension(self, dimkey, newdimvals, coordkey = None, **interpkwds):
        """
        Parameters
        ----------
        self : the file to interpolate from must have VGLVLS
        dimkey : the new dimension for interpolation
        newdimvals : the new values to interpolate to
        coordkey : the variable to use as the old coordinate values
        interptype : 'linear' or 'conserve'
             linear : uses a linear interpolation
             conserve : uses a mass conserving interpolation
        extrapolate : allow extrapolation beyond bounds with linear, default False
        fill_value : set fill value (e.g, nan) to prevent extrapolation or edge 
                     continuation
        
        Returns
        -------
        outf - ioapi_base PseudoNetCDFFile with al variables interpolated
        
        Notes
        -----
        When extrapolate is false, the edge values are used for points beyond
        the inputs.
        """
        from ..coordutil import getinterpweights

        if coordkey is None:
            olddimvals = self.variables[dimkey]
        else:
            olddimvals = self.variables[coordkey]
        if olddimvals.ndim == 1 and newdimvals.ndim == 1:
            weights = getinterpweights(olddimvals, newdimvals, **interpkwds)
            def interpd(data):
                if data.ndim == 1:
                    newdata = (weights * data[:, None]).sum(0)
                else:
                    newdata = (weights[None, :, :, None, None] * data[:, :, None]).sum(1)
                return newdata
            outf = self.applyAlongDimensions(**{dimkey: interpd})
        else:
            outf = self.copy(props = True, dimensions = False, variables = False)
            olddim = olddimvals.dimensions
            newdim = newdimvals.dimensions
            if olddim != newdim:
                raise ValueError('Can only interpolate if coordinate variable have the same named dimensions')
            dimaxis = olddim.index(dimkey)
            ndl = newdimvals.shape[dimaxis]
            for dk, dv in self.dimensions.items():
                if dk == dimkey:
                    dl = ndl
                else:
                    dl = len(dv)
                ndv = outf.createDimension(dk, dl)
                ndv.setunlimited(dv.isunlimited())
            
            for vk, vv in self.variables.items():
                nvv = outf.copyVariable(vv, key = vk, withdata = vv.dimensions != newdim)
            
            
            Ni, Nk = olddimvals.shape[:dimaxis], olddimvals.shape[dimaxis+1:]
            s_ = np.s_
            for ii in np.ndindex(Ni):
                for kk in np.ndindex(Nk):
                    od = olddimvals[ii + s_[:,] + kk]
                    nd = newdimvals[ii + s_[:,] + kk]
                    weights = getinterpweights(od, nd, **interpkwds)
                    for nvk, nvv in outf.variables.items():
                        if nvv.dimensions != newdim: continue
                        vv = self.variables[nvk]
                        nvv[ii + s_[...,] + kk] = (weights * vv[ii + s_[:,] + kk][:, None]).sum(0)
        return outf
            
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
                    if dvar.ndim != 1:
                        dvar = np.arange(len(dv))
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
    
    def getTimes(self, datetype = 'datetime', bounds = False):
        """
        Get an array of datetime objects
        
        Parameters
        ----------
        datetype : 'datetime' or numpy.dtype
        bounds : get time boundaries
        
        Returns
        -------
        array : array of datetime objects or array of numpy's datetype type
        Notes
        -----
        self must have a time or TFLAG variable
        """
        from PseudoNetCDF.coordutil import _parse_ref_date
        from datetime import date, datetime, timedelta, timezone
        utc = timezone.utc

        _calendaryearlike = {'noleap': 1970, '365_day': 1970, 'all_leap': 1972, '366_day': 1972}
               
        if 'time' in self.variables.keys():
            time = self.variables['time']
            timeunits = time.units.strip()
            calendar = getattr(time, 'calendar', 'gregorian').lower()
            if 'since' in timeunits:
                unit, base = timeunits.split(' since ')
                # Get the reference date
                refdate = _parse_ref_date(base)
                
                if calendar in _calendaryearlike:
                    refyear = refdate.year
                    # Get a year for relative day calculations
                    yearlike = _calendaryearlike[calendar]
                    # In that year, how many seconds and days are there
                    yearseconds = (date(yearlike + 1, 1, 1) - date(yearlike, 1, 1)).total_seconds()
                    yeardays = yearseconds / 3600 / 24
                    
                    # Get a new reference date in yearlike
                    crefdate = datetime(yearlike, 1 ,1, tzinfo = utc)
                    if refdate.month != 1 or refdate.day != 1:
                        # Get start date in yearlike
                        refcdate = datetime(yearlike, refdate.month, refdate.day, tzinfo = utc)
                        # Calculate delta in years
                        addyears = (crefdate - refcdate).total_seconds() / yearseconds
                    else:
                        addyears = 0
                    # Convert time to fractional years, including change in reference
                    fracyearincrs = time[:] / {'years': 1, 'days': yeardays, 'hours': yeardays*24, 'minutes': yeardays*24*60, 'seconds': yeardays*24*60}[unit] + addyears
                    # Split into years and days
                    yearincrs = np.array(fracyearincrs // 1).astype('i')
                    dayincrs = (fracyearincrs % 1) * yeardays
                    # Add days to the calendar year reference
                    cdays = [crefdate + timedelta(days = dayinc) for dayinc in dayincrs]
                    try:
                        # Combine calendar specific month and day with new year
                        out = np.array([datetime(refyear + yearinc, cday.month, cday.day, tzinfo = utc) for yearinc, cday in zip(yearincrs, cdays)])
                    except:
                        warn('Years calculated from %d day year, but month/days calculated for actual year. Usually means data has Feb 29th in a non leap year' % yeardays)
                        out = np.array([datetime(refyear + yearinc, 1, 1, tzinfo = utc) + timedelta(days = float(dayinc)) for yearinc, dayinc in zip(yearincrs, dayincrs)])

                else:
                    out = refdate + np.array([timedelta(**{unit: float(i)}) for i in time[:]])
                    
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
            out = np.array([datetime(yyyy, 1, 1, tzinfo = utc) + timedelta(days = day - 1) for yyyy, day in zip(yyyys, days)])
        elif hasattr(self, 'SDATE') and hasattr(self, 'STIME') and \
             hasattr(self, 'TSTEP') and 'TSTEP' in self.dimensions:
            refdate = datetime.strptime('%07d %06d+0000' % (self.SDATE, self.STIME), '%Y%j %H%M%S%z')
            tstepstr = '%06d' % self.TSTEP
            timeincr = timedelta(seconds = int(tstepstr[-2:])) + \
                       timedelta(minutes = int(tstepstr[-4:-2])) + \
                       timedelta(hours   = int(tstepstr[:-4]))
            ntimes = len(self.dimensions['TSTEP'])
            if bounds:
                ntimes += 1
            timeincrs = timeincr * np.arange(ntimes)
            out = refdate + timeincrs
        elif 'tau0' in self.variables.keys():
            out = datetime(1985, 1, 1, 0, tzinfo = utc) + np.array([timedelta(hours =i) for i in self.variables['tau0'][:]])
        else:
            raise ValueError('cannot understand time for file')
        if datetype == 'datetime':
            return out
        else:
            return np.array(out, dtype = datetype)

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
        isarray = {dk: not np.isscalar(dv) and not isinstance(dv, slice) for dk, dv in dimslices.items()}
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

    def __iter__(self):
        for k in self.keys():
            yield k
    
    def __len__(self):
        return len([k for k  in self.keys()])
    
    def items(self):
        return [(k, self[k]) for k in self.keys()]
    
    def __contains__(self, k):
        return k in [k for k in self.keys()]
    

if __name__ == '__main__':
    unittest.main()
