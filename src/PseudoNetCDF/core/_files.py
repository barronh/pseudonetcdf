__all__ = ['PseudoNetCDFFile', 'netcdf', 'PseudoNetCDFVariables']
import unittest
from PseudoNetCDF._getreader import registerreader
from PseudoNetCDF.netcdf import NetCDFFile, NetCDFVariable
from PseudoNetCDF.pncwarn import warn
from collections import OrderedDict
from ._dimensions import PseudoNetCDFDimension, PseudoNetCDFDimensions
from ._variables import PseudoNetCDFVariable, PseudoNetCDFMaskedVariable
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

    def __missing__(self, key):
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
        longname = '.'.join(
            [p for p in pieces[1:-1] if '_' != p[0] and p not in ('core',)] +
            [pieces[-1]])
        if len(cls.mro()) > 2:
            if name not in ('PseudoNetCDFFile', 'WrapPnc'):
                shortl = registerreader(name, cls)
                longl = registerreader(longname, cls)
                if not (shortl or longl):
                    warn('Not registered either as ' + name +
                         ' or ' + longname)
        super(PseudoNetCDFType, cls).__init__(name, bases, clsdict)


PseudoNetCDFSelfReg = PseudoNetCDFType('pnc', (object,), dict(__doc__='Test'))


class PseudoNetCDFFile(PseudoNetCDFSelfReg, object):
    """
    PseudoNetCDFFile provides an interface and standard set of
    methods that a file should present to act like a netCDF file
    using the Scientific.IO.NetCDF.NetCDFFile interface.
    """

    def getMap(self, maptype='basemap_auto', **kwds):
        """
        Description

        Parameters
        ----------
        maptype : string
            choices 'basemap', 'basemap_auto', 'cartopy' (not yet)
            basemap : attempts to open a basemap with only supplied kwds
            basemap_auto : automatically adds llcrnrlon,llcrnrlat,u
                           rcrnrlon,urcrnrlat based on longitude_bounds
        **kwds : keywords
            for basemap or cartopy

        Returns
        -------
        map : basemap or cartopy axis
        """
        if maptype.startswith('basemap'):
            from PseudoNetCDF.coordutil import basemap_from_proj4
            if (
                maptype.endswith('_auto') and
                'longitude_bounds' in self.variables and
                'latitude_bounds' in self.variables
            ):
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
                edges = dict(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                             urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
                kwds.update(edges)
            myproj = self.getproj(withgrid=True, projformat='proj4')
            return basemap_from_proj4(myproj, **kwds)
        elif maptype == 'cartopy':
            raise ValueError('cartopy is not yet implemented')
        else:
            raise ValueError(
                'maptype must be basemap, basemap_auto, or cartopy')

    def getproj(self, withgrid=False, projformat='pyproj'):
        """
        Description

        Parameters
        ----------
        withgrid : boolean
            use grid units instead of meters
        projformat : string
            'pyproj' (default), 'proj4' or 'wkt' allows function to
            return a pyproj projection object or a string in the
            format of proj4 or WKT

        Returns
        -------
        proj : string pyproj.Proj
             (wkt, proj4) or pyprojProj (pyproj)
        """
        if projformat == 'pyproj':
            from PseudoNetCDF.coordutil import getproj
            return getproj(self, withgrid=withgrid)
        elif projformat == 'proj4':
            from PseudoNetCDF.coordutil import getproj4
            return getproj4(self, withgrid=withgrid)
        elif projformat == 'wkt':
            from PseudoNetCDF.coordutil import getprojwkt
            return getprojwkt(self, withgrid=withgrid)
        else:
            raise ValueError('projformat must be pyproj, proj4 or wkt')

    def ll2xy(self, lon, lat):
        """
        Converts lon/lat to x distances (no false easting/northing)

        Parameters
        ----------
        lon : scalar or iterable
            longitudes in decimal degrees
        lat : scalar or iterable
            latitudes in decimal degrees

        Returns
        -------
        x, y : tuple of arrays
            coordinates in map projection (meters or radians)
        """
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        return self.getproj()(lon, lat)

    def _getzdim(self):
        """
        Heuristic function finds and returns the name of the vertical dimension
        based on common dimension names. Overwrite method for explicit control.

        Returns
        -------
        result : {'LAY', 'bottom_top', 'level', 'lev', 'layer', 'lay', 'z'}

        Examples
        --------
        >>> f = pnc.pncopen('', format='PseudoNetCDFFile')
        >>> f.createDimension('t', 4)
        >>> f.createDimension('x', 4)
        >>> f.createDimension('y', 4)
        >>> f.createDimension('z', 4)
        >>> f._getzdim()
        'z'
        """
        for dk in 'LAY bottom_top level lev layer lay z'.split():
            if dk in self.dimensions:
                return dk
        else:
            raise KeyError('Could not find z dimensions')

    def _gettdim(self):
        """
        Heuristic function finds and returns the name of the time dimension
        based on common dimension names. Overwrite method for explicit control.

        Returns
        -------
        result : {'time', 'TSTEP', 't', 'Time'}

        Examples
        --------
        >>> f = pnc.pncopen('', format='PseudoNetCDFFile')
        >>> f.createDimension('t', 4)
        >>> f.createDimension('x', 4)
        >>> f.createDimension('y', 4)
        >>> f.createDimension('z', 4)
        >>> f._gettdim()
        't'
        """
        for dk in 'time TSTEP t Time'.split():
            if dk in self.dimensions:
                return dk
        else:
            raise KeyError('Could not find t dimensions')

    def _getydim(self):
        """
        Heuristic function finds and returns the name of the y or width
        dimension based on common dimension names. Overwrite method for
        explicit control.

        Returns
        -------
        result : {'latitude', 'lat', 'south_north', 'ROW', 'y'}

        Examples
        --------
        >>> f = pnc.pncopen('', format='PseudoNetCDFFile')
        >>> f.createDimension('t', 4)
        >>> f.createDimension('x', 4)
        >>> f.createDimension('y', 4)
        >>> f.createDimension('z', 4)
        >>> f._getydim()
        'y'
        """
        for dk in 'latitude lat south_north ROW y'.split():
            if dk in self.dimensions:
                return dk
        else:
            raise KeyError('Could not find y dimensions')

    def _getxdim(self):
        """
        Heuristic function finds and returns the name of the x or height
        dimension based on common dimension names. Overwrite method for
        explicit control.

        Returns
        -------
        result : {'longitude', 'lon', 'west_east', 'COL', 'x'}

        Examples
        --------
        >>> f = pnc.pncopen('', format='PseudoNetCDFFile')
        >>> f.createDimension('t', 4)
        >>> f.createDimension('x', 4)
        >>> f.createDimension('y', 4)
        >>> f.createDimension('z', 4)
        >>> f._getydim()
        'y'
        """
        for dk in 'longitude long lon west_east COL x'.split():
            if dk in self.dimensions:
                return dk
        else:
            raise KeyError('Could not find x dimensions')

    def ll2ij(self, lon, lat, bounds='ignore', clean='none'):
        """
        Converts lon/lat to 0-based indicies (0,M), (0,N)

        Parameters
        ----------
        lon : scalar or iterable
            longitudes in decimal degrees
        lat : scalar or iterable
            latitudes in decimal degrees
        bounds : string
            ignore, error, warn if i,j are out of domain
        clean : string
            none - return values regardless of bounds;
            mask - mask values out of bounds;
            clip - return min(max(0, v), nx - 1)

        Returns
        -------
        i, j : indices (0-based) for variables
        """
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        p = self.getproj(withgrid=True)
        x, y = p(lon, lat)
        i = np.asarray(x).astype('i')
        j = np.asarray(y).astype('i')
        nx = len(self.dimensions[self._getxdim()])
        ny = len(self.dimensions[self._getydim()])
        if bounds == 'ignore':
            pass
        else:
            lowi = (i < 0)
            lowj = (j < 0)
            highi = (i >= nx)
            highj = (j >= ny)
            outb = (lowi | lowj | highi | highj)
            nout = outb.sum()
            if nout > 0:
                message = '{} Points out of bounds; {}'.format(
                    nout, np.where(outb))
                if bounds == 'error':
                    raise ValueError(message)
                else:
                    warn(message)
        if clean == 'clip':
            i = np.minimum(np.maximum(i, 0), nx - 1)
            j = np.minimum(np.maximum(j, 0), ny - 1)
        elif clean == 'mask':
            i = np.ma.masked_greater(np.ma.masked_less(i, 0), nx - 1)
            j = np.ma.masked_greater(np.ma.masked_less(j, 0), ny - 1)

        return i, j

    def xy2ll(self, x, y):
        """
        Converts x, y to lon, lat (no false easting/northing)

        Parameters
        ----------
        x : scalar or iterable
            projected west-east coordinates
        y : scalar or iterable
            projected south-north coordinates

        Returns
        -------
        lon, lat : scalars or iterables
            longitudes and latitudes in decimal degrees
        """
        x = np.asarray(x)
        y = np.asarray(y)
        p = self.getproj()
        lon, lat = p(x, y, inverse=True)
        return lon, lat

    def ij2ll(self, i, j):
        """
        Converts i, j to lon, lat (no false easting/northing)
        using cell centers assuming 0-based i/j

        Parameters
        ----------
        i : scalar/iterable
            indicies (0-based) for the west-east dimension
        j : scalar/iterable
            indicies (0-based) for the south-north dimension

        Returns
        -------
        lon, lat : scalars or iterables
            longitudes and latitudes in decimal degrees
        """
        i = np.asarray(i)
        j = np.asarray(j)
        p = self.getproj(withgrid=True)
        lon, lat = p(i + 0.5, j + 0.5, inverse=True)
        return lon, lat

    def date2num(self, time, timekey='time'):
        """
        Parameters
        ----------
        time : array-like
            array of datetime.datetime objects
        timekey : str
            time variable key which requires units and should have calendar.
            If calendar is missing, standard is the default. default 'time'

        Returns
        -------
        num : array-like
            time in relative time as defined by units of time variable
            (i.e., timekey) which defaults to 'time'
        """
        from netCDF4 import date2num
        try:
            from datetime import timezone
            utc = timezone.utc
        except ImportError:
            from datetime import tzinfo
            utc = tzinfo.utc

        time = np.asarray(time)
        # netCDF4 date2num is timezone naive; assumes UTC when not
        # specified and converts to UTC internally
        # so, if a tzinfo is involved, it should be removed
        if any([t.tzinfo is not None for t in time[:]]):
            time = np.array([
                t.astimezone(utc).replace(tzinfo=None) for t in time[:]
            ])

        timeunits = self.variables[timekey].units.strip()
        calendar = getattr(self.variables[timekey], 'calendar', 'standard')
        num = date2num(time, timeunits, calendar.strip())
        return num

    def time2idx(self, time, dim='time', timekey=None, **kwds):
        """
        Convert datetime objects to dimension indices

        Parameters
        ----------
        time : array-like
            array of datetime.datetime objects
        dim : str
            dimension name for val2idx
        timekey : str
            time variable key. None defaults to dim
        kwds : mappable
            see val2idx

        Returns
        -------
        idx : array-like
            time index (0-based)
        """
        if timekey is None:
            timekey = dim

        time = np.asarray(time)
        nums = self.date2num(time, timekey=timekey)
        return self.val2idx(dim=dim, val=nums, **kwds)

    def time2t(self, time, ttype='nearest', index=True):
        """
        Parameters
        ----------
        time : array of datetime.datetime objects
        interp : 'nearest', 'bounds', 'bounds_close'
        index : return index

        Returns
        -------
        t : fractional time or if index, integers for indexing
        """
        warn(
            "time2t is deprecated; transition to time2idx or date2num",
            DeprecationWarning
        )
        time = np.asarray(time)
        if ttype not in ('nearest', 'bounds', 'bounds_close'):
            warn('{} is not an option; defaulting to nearest'.format(ttype))
            ttype = 'nearest'
        if ttype == 'nearest':
            mytimes = self.getTimes()
        else:
            mytimes = self.getTimes(bounds=True)
        idx = np.arange(mytimes.size)

        # Find minimum resoultion
        for res in [
            'microsecond', 'second', 'minute', 'hour', 'day', 'month', 'year'
        ]:
            minres = res
            resvals = [getattr(d, res) for d in time]
            myresvals = [getattr(d, res) for d in mytimes]
            if not np.allclose(myresvals, 0):
                break
            if not np.allclose(resvals, 0):
                break

        # Translate minimum resolution to numpy datetime
        if minres == 'microsecond':
            tu = 'datetime64[ns]'
        elif minres == 'second':
            tu = 'datetime64[s]'
        elif minres == 'minute':
            tu = 'datetime64[m]'
        elif minres == 'hour':
            tu = 'datetime64[h]'
        elif minres in ('year', 'month', 'day'):
            tu = 'datetime64[D]'
        else:
            tu = 'datetime64[ns]'

        # Convert input time to numpy datetime at resolution
        x = time.astype(tu).astype('d')
        # Convert file's time to numpy datetime at resolution
        xp = mytimes.astype(tu).astype('d')

        # Use interpolation methods with no bounding for nearest
        # and bounds_close
        if ttype in ('nearest', 'bounds_close'):
            out = np.interp(x, xp, idx)
            if index:
                imin = 0
                imax = idx[-1] + (0 if ttype == 'nearest' else -1)
                out = np.minimum(np.maximum(out, imin), imax)
                if ttype == 'nearest':
                    out = np.round(out, 0).astype('i')
                else:
                    out = np.floor(out).astype('i')
        # Use interpolation methods with bounding for nearest
        else:
            out = np.interp(x, xp, idx, left=np.nan, right=np.nan)
            if index:
                out = np.ma.masked_less(np.ma.floor(out).astype('i'), 0)

        return out

    def val2idx(
        self, dim, val,
        method='nearest', bounds='warn', left=None, right=None, clean='mask'
    ):
        """
        Convert coordinate values to indices

        Parameters
        ----------
        dim : str
            name of dimensions, which must have a coordinate variable
        val : array-like
            value in coordinate space
        method : str
            nearest, bounds, exact - each calculates the index differently
             - nearest : uses interp with coord values and rounds
             - bounds : uses interp between bounding values and truncates
             - exact : returns indices for exact coord values with other
                       indices masked (clean keyword has no effect)
        bounds : str
            ignore, error, warn if i,j are out of domain
        left : scalar
            see np.interp
        right : scalar
            see np.interp
        clean : {'none', 'mask'}
            none - return values regardless of bounds;
            mask - mask invalid values (use with left/right=np.nan);
                   has no affect with method exact

        Returns
        -------
        i : array-like
            indices (0-based) for variables
        """
        val = np.asarray(val)
        if method not in ('exact', 'nearest', 'bounds'):
            raise NotImplementedError(method + ' method not implemented')

        if bounds not in ('ignore', 'warn', 'error'):
            raise NotImplementedError(bounds + ' bounds not implemented')

        if clean not in ('none', 'mask'):
            raise NotImplementedError(clean + ' clean not implemented')

        dimv = self.variables[dim]
        dimvals = dimv[...]
        if dimvals.ndim > 1:
            raise ValueError(
                'val2idx is only implemented for 1-D coordinate variables'
            )
        if method == 'bounds':
            bounds_keys = [dim + '_bounds', dim + '_bnds']
            if hasattr(dimv, 'bounds'):
                bounds_keys.insert(0, dimv.bounds)

            for dimbk in bounds_keys:
                if dimbk not in self.variables:
                    continue
                dimbv = self.variables[dimbk]
                if dimbv.ndim == 1:
                    dimevals = dimbv[:]
                elif dimbv.ndim == 2 and dimbv.shape[1] == 2:
                    dimevals = np.append(dimbv[:, 0], dimbv[-1, 1])
                else:
                    raise ValueError(
                        'val2idx is only implemented for 1-D or 2-D bounds' +
                        '; {} has {}'.format(dimbk, dimbv.shape)
                    )
                break
            else:
                warn('Approximating bounds for val2idx {}'.format(dim))
                dval = np.diff(dimvals) / 2
                start = dimvals[:1]
                end = dimvals[-1:]
                if (dval == dval[0]).all():
                    start -= dval[0]
                    end += dval[-1]

                dimevals = np.concatenate([
                    start,
                    dimvals[1:] - dval,
                    end
                ])
        else:
            dimevals = dimvals

        idx = np.arange(dimevals.size)
        ddimevals = np.diff(dimevals)

        if (ddimevals < 0).all():
            dimevals[::-1]
            idx = idx[::-1]
        elif (ddimevals > 0).all():
            pass
        else:
            raise ValueError('coordinate is neither ascending nor descending')

        # import pdb; pdb.set_trace()
        # left = dimevals[0] - 1
        # right = dimevals[-1] + 1
        fidx = np.interp(val, dimevals, idx, left=left, right=right)
        if method == 'bounds':
            if right is None or right == dimevals[-1]:
                fidx = np.minimum(fidx, dimvals.size - 1)

        if method == 'exact':
            fidx = np.ma.masked_where(~np.in1d(val, dimvals), fidx)

        if clean == 'mask':
            outfidx = np.ma.masked_invalid(fidx)
        else:
            outfidx = fidx

        if bounds != 'ignore':
            isleft = val < dimevals[0]
            isright = val > dimevals[-1]
            isout = isleft | isright
            outmesg = 'Values are out of bounds:\n{}'.format(val[isout])
            if bounds == 'warn':
                warn(outmesg)
            elif bounds == 'error':
                raise ValueError(bounds)

        if method == 'nearest':
            outidx = np.round(outfidx, 0).astype('i')
        else:
            outidx = outfidx.astype('i')

        return outidx

    def get_dest(self):
        """
        Returns
        -------
        path : str
            path where a new file is created on some action

        Notes
        -----
        If None, a file is created in memory.
        Else, a netcdf file is created on disk.
        """
        return getattr(self, '_destination', None)

    def set_dest(self, path, **options):
        """
        Sets the path where a new file is created on some action

        Arguments
        ---------
        path : str
            path for new file
        **options : keywords for constructor
            options for new file creation

        Returns
        -------
        None
        """
        options['filename'] = path
        return object.__setattr__(self, '_destination', options)

    def get_varopt(self):
        """
        Get options

        Arguments
        ---------
        None

        Returns
        -------
        options : dictionary of options
        """
        return getattr(self, '_varopt', {}).copy()

    def set_varopt(self, **options):
        """
        Set options to be used when creating any Variable

        Arguments
        ---------
        **options : options for createVariable
            optional keywords to be supplied when creating new variables in
            destination file

        Returns
        -------
        None
        """
        return object.__setattr__(self, '_varopt', options)

    def _newlike(self):
        """
        Internal function to return a file of the same class if a
        PseudoNetCDFFile
        """
        if self.get_dest() is not None:
            outf = netcdf(**self.get_dest())
        elif isinstance(self, PseudoNetCDFFile):
            outt = type(self)
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        outf.set_varopt(**self.get_varopt())
        return outf

    def renameVariable(self, oldkey, newkey, inplace=False, copyall=True):
        """
        Rename variable (oldkey)

        Parameters
        ----------
        oldkey : string
            variable to be renamed
        newkey : string
            new dame for variable
        inplace : boolean
            create the new variable in this netcdf file (default False)
        copyall : boolean
            if not inplace, should all variables be copied to new file

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with renamed variable (this file if inplace = True)
        """
        return self.renameVariables(**{oldkey: newkey})

    def renameVariables(self, inplace=False, copyall=True, **newkeys):
        """
        Rename variables for each oldkey: newkey dictionary item

        Parameters
        ----------
        **newkeys : dictionary
            where key is the oldkey and value is the newkey
        inplace : boolean
            create the new variable in this netcdf file (default False)
        copyall :boolean
            if not inplace, should all variables be copied to new file

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with renamed variable (this file if inplace = True)
        """
        if inplace:
            outf = self
        else:
            outf = self._copywith(
                props=True, dimensions=True, variables=copyall, data=copyall)

        for oldkey, newkey in newkeys.items():
            outf.copyVariable(self.variables[oldkey], key=newkey)
            if oldkey in outf.variables:
                del outf.variables[oldkey]

        return outf

    def renameDimension(self, oldkey, newkey, inplace=False):
        """
        Rename dimension (oldkey) in dimensions and in all variables

        Parameters
        ----------
        oldkey : string
            dimension to be renamed
        newkey : string
            new dame for dimension
        inplace : boolean
            create the new variable in this netcdf file (default False)

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with renamed variable (this file if inplace = True)
        """
        return self.renameDimensions(inplace=inplace, **{oldkey: newkey})

    def renameDimensions(self, inplace=False, **newkeys):
        """
        Rename dimension (oldkey) in dimensions and in all variables

        Parameters
        ----------
        **newkeys : dictionary
            where key is the oldkey and value is the newkey
        inplace : boolean
            create the new variable in this netcdf file (default False)

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with renamed variable (this file if inplace = True)
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

    def insertDimension(self, newonly=True, multionly=False, before=None,
                        after=None, inplace=False, **newdims):
        """
        Insert dimensions with keys and lengths from newdims


        Parameters
        ----------
        **newdims : dictionary
            where key is the new dimension and value is the length
        newonly : boolean
            Only add dimension to variables that do not already have it,
            default True
        multionly : boolean
            Only add dimension if there are already more than one (good
            for ignoring coordinate dimensions)
        before : string
            if variable has this dimension, insert the new dimension before
            it. Otherwise, add to the beginning. (before takes precedence)

        after : string
            if variable has this dimension, insert the new dimension after it.
            Otherwise, add to the beginning.

        inplace : boolean
            create the new variable in this netcdf file (default False)

        Returns
        -------
        outf : PseudoNetCDFFile
            instance will new dimension in dimensions and variables

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
            outf = self.copy(variables=False)

        for dk, dv in newdims.items():
            if dk not in outf.dimensions:
                outf.createDimension(dk, dv)

            for vk, vv in self.variables.items():
                vdims = list(vv.dimensions)
                if (
                    (newonly and dk in vdims) or
                    (multionly and len(vdims) == 1)
                ):
                    outf.copyVariable(vv, key=vk, withdata=True)
                    continue
                ndims = [_dk for _dk in vdims]
                if before in vdims:
                    bi = vdims.index(before)
                elif after in vdims:
                    bi = vdims.index(after) + 1
                elif not (before is None and after is None):
                    outf.copyVariable(vv, key=vk, withdata=True)
                    continue
                else:
                    bi = 0
                ndims.insert(bi, dk)
                var = outf.copyVariable(
                    vv, key=vk, dimensions=ndims, withdata=False)

                var[...] = np.expand_dims(vv[...], axis=bi)
        return outf

    def reorderDimensions(self, oldorder, neworder, inplace=False):
        """
        Evaluate expr and return a PseudoNetCDFFile object with resutl

        Parameters
        ----------
        oldorder : iterable of strings
            dimension names in existing order
        neworder : iterable of strings
            dimension names in new order

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with dimensions reordered in variables
        """
        if inplace:
            outf = self
        else:
            outf = self.copy(variables=True)
        oldorder = tuple(oldorder)
        neworder = tuple(neworder)
        for vk, vv in self.variables.items():
            varneworder = [dk for dk in neworder if dk in vv.dimensions]
            varorder = [dk for dk in vv.dimensions]
            if len(varneworder) > 0:
                newvals = vv[:].copy()
                for newdi, newdk in enumerate(varneworder):
                    axisidx = varorder.index(newdk)
                    if axisidx == newdi:
                        continue
                    newvals = np.rollaxis(newvals, axis=axisidx, start=newdi)
                    varorder.pop(axisidx)
                    varorder.insert(newdi, newdk)
                assert(varorder == varneworder)
                newvals.dimensions = tuple(varorder)
                outf.variables[vk] = newvals
            else:
                pass

        return outf

    def eval(self, expr, inplace=False, copyall=False):
        """
        Evaluate expr and return a PseudoNetCDFFile object with resutl

        Parameters
        ----------
        expr : string
            expression to evaluate
        inplace : boolean
            create the new variable in this netcdf file (default False)
        copyall : boolean
            if not inplace, should all variables be copied to new file

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with renamed variable (this file if inplace = True)
        """
        import numpy as np
        # Copy file to temporary PseudoNetCDF file
        comp = compile(expr, 'none', 'exec')
        vardict = {k: v for k, v in self.variables.items()}
        for pk in self.ncattrs():
            if pk not in vardict:
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
            if copyall:
                outf = self.copy()
            else:
                newkeys = [key]
                outf = self.subsetVariables(newkeys)
                try:
                    del outf.variables[key]
                except Exception:
                    pass

        propd = dict([(k, getattr(tmpvar, k)) for k in tmpvar.ncattrs()])
        propd['expression'] = expr
        dimt = tmpvar.dimensions
        vardict['outf'] = self

        # Assign expression to new variable.
        exec(comp, None, vardict)
        assignedkeys = [s.get_name() for s in symbols if s.is_assigned()]
        assignedkeys = [k for k in assignedkeys if k in vardict]
        for key in assignedkeys:
            try:
                del outf.variables[key]
            except Exception:
                pass
            val = vardict[key]
            # if the output variable has no dimensions,
            # there is likely a problem and the output
            # should be defined.
            if (
                isinstance(val, (PseudoNetCDFVariable,)) and
                val.dimensions != ()
            ):
                outf.variables[key] = val
            else:
                outf.createVariable(key, val.dtype.char,
                                    dimt, values=val, **propd)

        return outf

    def mask(
        self, where=None, less=None, less_equal=None, greater=None,
        greater_equal=None, values=None, equal=None, invalid=False, mask=None,
        dims=None, fill_value=-999, coords=False, verbose=0
    ):
        """
        Apply mask to all variables of same shape or just where dimensions
        match.

        Parameters
        ----------
        where : array-like
            boolean array to use as a mask see numpy.ma.masked_where
        greater : scalar
            mask when values are greater than this value see
            numpy.ma.masked_greater
        less : scalar
            mask when values are less than this value see
            numpy.ma.masked_less
        greater_equal : scalar
            mask when values are greater than or equal to this value see
            numpy.ma.masked_greater
        less_equal : scalar
            mask when values are less than or equal to this value see
            numpy.ma.masked_less
        values : scalar
            mask when values are equal to this value within standard floating
            point see numpy.ma.masked_values
        equal : scalar
            mask when values are exactly this value (i.e., integers) see
            numpy.ma.masked_equal
        invalid : boolean
            mask when values are invalid, see numpy.ma.masked_invalid
        mask : array-like
            alias for where
        dims : iterable of strings
            only apply "mask" or "where" to variables with these dimensions
            no effet on other masks
        fill_value : scalar
            value to use as the fill_value for new arrays
        coords : boolean
            if True, apply masks to coordinate variables. Default, False
        verbose : int
            level of verbosity for function, mostly for debugging

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with masked variables

        See Also
        --------
        numpy.ma : all masks are passing throught to numpy.ma.masked_...

        Notes
        -----
        mask options are not mutually exclusive. the order is where, greater,
        greater_equal, less, less_equal, values, equal, invalid
        """
        from time import time
        if where is None and mask is not None:
            where = mask
        elif mask is not None and where is not None:
            raise ValueError(
                'mask is an alias to where; supply one or the other, ' +
                'but not both'
            )

        if dims is None:
            maskdims = getattr(where, 'dimensions', dims)
        else:
            maskdims = dims

        coordkeys = self.getCoords()
        outf = self.copy(variables=False)
        nitems = len(self.variables.keys())
        if verbose == 1:
            print('|' + '=' * nitems + '|', flush=True)
            print('|', end='', flush=True)
        for vk, vv in self.variables.items():
            if verbose > 1:
                t0 = time()
                print(vk, end='', flush=True)
            elif verbose > 0:
                print('.', end='', flush=True)
            newvar = outf.copyVariable(
                vv, key=vk, fill_value=fill_value, withdata=False
            )
            if vk in coordkeys and not coords:
                newvar[...] = vv[...]
                continue

            vals = vv[...]
            if where is not None:
                if (
                    maskdims == vv.dimensions or
                    (
                        maskdims is None and
                        where.shape == vals.shape
                    )
                ):
                    vals = np.ma.masked_where(where, vals)
                    # np.ma.putmask(newvar[...], mask==False, vv[...])
                    # vals = newvar[...]
                # else:
                #     vals = vv[...]

            if greater is not None:
                vals = np.ma.masked_greater(vals, greater)

            if greater_equal is not None:
                vals = np.ma.masked_greater_equal(vals, greater_equal)

            if less is not None:
                vals = np.ma.masked_less(vals, less)

            if less_equal is not None:
                vals = np.ma.masked_less_equal(vals, less_equal)

            if values is not None:
                vals = np.ma.masked_values(vals, values)

            if equal is not None:
                vals = np.ma.masked_equal(vals, equal)

            if invalid:
                vals = np.ma.masked_invalid(vals)

            if verbose > 1:
                t1 = time()
                print(t1 - t0, end='\n', flush=True)

            newvar[...] = vals[...]
        if verbose > 1:
            pass
        elif verbose > 0:
            print('')
        return outf

    def plot(self, varkey, plottype=None, ax_kw=None,
             plot_kw=None, cbar_kw=None, map_kw=None, dimreduction='mean'):
        """
        Parameters
        ----------
        varkey : string
            the variable to plot
        plottype : string
            any dimension name pair delimited by a hyphen (e.g.,
            longitude-latitude, latitude-pressure, longitude-pressure,
            vertical-profile, time-longitude, time-latitude, time-pressure)
            defaults to the last two dimensions.
        ax_kw : dictionary
            keywords for the axes to be created
        plot_kw : dictionary
            keywords for the plot (plot, scatter, or pcolormesh) to be
            created
        cbar_kw : dictionary or bool or None
            keywords for the colorbar; if True or None, use defaults.
            If False, do not create a colorbar
        map_kw : dictionary or bool or None
            keywords for the getMap routine, which is only used with
            map capable dimensions (ie, plottype='longitude-latitude')
            If True or None, use default configuration dict(countries=True,
            coastlines=True, states=False, counties=False). If False,
            do not draw a map.
        dimreduction : string or function
            dimensions not being used in the plot are removed using
            applyAlongDimensions(dimkey=dimreduction) where each dimenions
        """
        import matplotlib.pyplot as plt
        from ..coordutil import getbounds

        if ax_kw is None:
            ax_kw = {}

        if plot_kw is None:
            plot_kw = {}

        if cbar_kw is None or cbar_kw is True:
            cbar_kw = {}

        if map_kw is None or map_kw is True:
            map_kw = {}
        elif map_kw is not False:
            map_kw = map_kw.copy()

        apply2dim = {}
        var = self.variables[varkey]
        if plottype is None:
            vdims = var.dimensions
            if len(vdims) > 1:
                plottype = '-'.join(var.dimensions[-2:][::-1])
            else:
                plottype = var.dimensions[0] + '-profile'

        varunit = varkey
        if hasattr(var, 'units'):
            varunit += ' ' + var.units.strip()

        dimlens = dict([(dk, len(self.dimensions[dk]))
                        for dk in var.dimensions])
        dimpos = dict([(dk, di) for di, dk in enumerate(var.dimensions)])
        xkey, ykey = plottype.split('-')
        if not ykey == 'profile':
            for dimkey in list(dimlens):
                if dimkey not in (xkey, ykey) and dimlens.get(dimkey, 1) > 1:
                    apply2dim[dimkey] = dimreduction

        if len(apply2dim) > 0:
            subsetkeys = [varkey]
            if xkey in self.variables:
                subsetkeys.append(xkey)
            if ykey in self.variables:
                subsetkeys.append(ykey)
            myf = self.subsetVariables(subsetkeys)\
                      .applyAlongDimensions(**apply2dim)
            var = myf.variables[varkey]
            dimlens = dict([(dk, len(self.dimensions[dk]))
                            for dk in var.dimensions])
        else:
            myf = self
        if ykey in ('profile',):
            vaxi = var.dimensions.index(xkey)
            vsize = var.shape[vaxi]
            vals = np.rollaxis(var[:], vaxi).reshape(vsize, -1)
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
            x1 = vals[:].min(1)
            xm = vals[:].mean(1)
            x2 = vals[:].max(1)
            ax.fill_betweenx(y=y, x1=x1, x2=x2, label=varkey + '(min, max)')
            ax.plot(xm, y, label=varkey, **plot_kw)
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
        if cbar_kw is not False:
            cbar_kw.setdefault('label', varunit)
            ax.figure.colorbar(p, **cbar_kw)
        if xkey == 'time':
            ax.xaxis.set_major_formatter(
                plt.matplotlib.dates.AutoDateFormatter(
                    plt.matplotlib.dates.AutoDateLocator()
                )
            )
        if ykey == 'time':
            ax.yaxis.set_major_formatter(
                plt.matplotlib.dates.AutoDateFormatter(
                    plt.matplotlib.dates.AutoDateLocator()
                )
            )
        mappabledims = (
            plottype == 'longitude-latitude' or
            plottype == 'lon-lat'
        )
        if mappabledims and map_kw is not False:
            try:
                map_kw = map_kw.copy()
                coastlines = map_kw.pop('coastlines', True)
                countries = map_kw.pop('countries', True)
                states = map_kw.pop('states', False)
                counties = map_kw.pop('counties', False)
                bmap = myf.getMap(**map_kw)
                if coastlines:
                    bmap.drawcoastlines(ax=ax)
                if countries:
                    bmap.drawcountries(ax=ax)
                if states:
                    bmap.drawstates(ax=ax)
                if counties:
                    bmap.drawcounties(ax=ax)
            except Exception:
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
        attdict : dictionary
            key/value pairs of properties

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
        attdict : dictionary
            key/value pairs of properties
        """
        outd = OrderedDict()
        for pk in self.ncattrs():
            outd[pk] = getattr(self, pk)
        return outd

    @classmethod
    def from_ncvs(cls, *invars, **invarkw):
        """
        Arguments
        ---------
        invars : list
            NetCDF-like variable must have standard_name, long_name or name

        invarkw : kwds
            NetCDF-like variables

        Returns
        -------
        outf : PseudoNetcdf-like file
        """
        outf = cls()
        for invar in invars:
            for dk, ds in zip(invar.dimensions, invar.shape):
                if dk in outf.dimensions:
                    assert(ds == len(outf.dimensions[dk]))
                else:
                    outf.createDimension(dk, ds)
            outf.copyVariable(invar)

        for inkey, invar in invarkw.items():
            for dk, ds in zip(invar.dimensions, invar.shape):
                if dk in outf.dimensions:
                    assert(ds == len(outf.dimensions[dk]))
                else:
                    outf.createDimension(dk, ds)
            outf.copyVariable(invar, key=inkey)

        return outf

    @classmethod
    def from_ncf(cls, infile):
        """
        Arguments
        ---------
        infile : PseudoNetCDF-like file

        Returns
        -------
        outf : PseudoNetcdf-like file
        """
        outf = cls()
        for pk in infile.ncattrs():
            pv = getattr(infile, pk)
            setattr(outf, pk, pv)

        for dk, dv in infile.dimensions.items():
            outf.copyDimension(dv, key=dk)

        for vk, vv in infile.variables.items():
            outf.copyVariable(vv, key=vk)

        return outf

    def _copywith(self, props=True, dimensions=True, variables=False,
                  data=False):
        """
        Internal function for making copies of the same type

        Parameters
        ----------
        props : boolean
            include properties (default: True)
        dimensions : boolean
            include dimensions (default: True)
        variables : boolean
            include variable structures (default: False)
        data : boolean
            include variable data (default: False)

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with some parts

        Notes
        -----
        Internal function does not return variables by default.
        This is useful for functions like slice, apply, eval, etc.

        The _ in _copywith means this is a private function and the
        call interface may change.
        """
        outf = self._newlike()
        outf._operator_exclude_vars = tuple(self._operator_exclude_vars)
        if props:
            for pk in self.ncattrs():
                setattr(outf, pk, getattr(self, pk))
        if dimensions:
            for dk, dv in self.dimensions.items():
                outf.copyDimension(dv, key=dk)

        if variables:
            for vk, vv in self.variables.items():
                outf.copyVariable(vv, key=vk, withdata=data)

        return outf

    def copy(self, props=True, dimensions=True, variables=True, data=True):
        """
        Function for making copies of the same type

        Parameters
        ----------
        props : boolean
            include properties (default: True)
        dimensions : boolean
            include dimensions (default: True)
        variables : boolean
            include variable structures (default: True)
        data : boolean
            include variable data (default: True)

        Returns
        -------
        outf : PseudoNetCDFFile instance
        """
        return self._copywith(props=props, dimensions=dimensions,
                              variables=variables, data=data)

    def interpDimension(self, dimkey, newdimvals, coordkey=None, **interpkwds):
        """
        Parameters
        ----------
        dimkey : string
            the new dimension for interpolation
        newdimvals : iterable
            the new values to interpolate to
        coordkey : string
            the variable to use as the old coordinate values
        interptype : string
            'linear' or 'conserve'. linear uses a linear interpolation
             conserve uses a mass conserving interpolation
        extrapolate : boolean
            allow extrapolation beyond bounds with linear, default False
        fill_value : numeric value
            set fill value (e.g, nan) to prevent extrapolation or edge
            continuation

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with all variables interpolated

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
                    newdata = (weights[None, :, :, None, None] *
                               data[:, :, None]).sum(1)
                return newdata
            outf = self.applyAlongDimensions(**{dimkey: interpd})
        else:
            outf = self.copy(props=True, dimensions=False, variables=False)
            olddim = olddimvals.dimensions
            newdim = newdimvals.dimensions
            if olddim != newdim:
                raise ValueError('Can only interpolate if coordinate ' +
                                 'variable have the same named dimensions')
            dimaxis = olddim.index(dimkey)
            ndl = newdimvals.shape[dimaxis]
            for dk, dv in self.dimensions.items():
                if dk == dimkey:
                    dl = ndl
                else:
                    dl = len(dv)
                outf.copyDimension(dv, key=dk, dimlen=dl)

            for vk, vv in self.variables.items():
                nvv = outf.copyVariable(
                    vv, key=vk, withdata=vv.dimensions != newdim)

            Ni, Nk = olddimvals.shape[:dimaxis], olddimvals.shape[dimaxis + 1:]
            s_ = np.s_
            for ii in np.ndindex(Ni):
                for kk in np.ndindex(Nk):
                    od = olddimvals[ii + s_[:, ] + kk]
                    nd = newdimvals[ii + s_[:, ] + kk]
                    weights = getinterpweights(od, nd, **interpkwds)
                    for nvk, nvv in outf.variables.items():
                        if nvv.dimensions != newdim:
                            continue
                        vv = self.variables[nvk]
                        interpedv = (weights *
                                     vv[ii + s_[:, ] + kk][:, None]).sum(0)
                        nvv[ii + s_[..., ] + kk] = interpedv
        return outf

    def applyAlongDimensions(self, verbose=0, **dimfuncs):
        """
        Similar to numpy.apply_along_axis, but for damed dimensions and
        processes dimensions as well as variables

        Parameters
        ----------
        dimfuncs : dictionary
            key value pairs where the key is a dimensions and the value is a
            1D function (func1d) or a dictionary. If the value is a dictionary
            it must include func1d as a function and any keyword arguments as
            additional options
        verbose : integer
            0 silent, 1 show variable, 2 show dimensions and variables

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with variables and dimensions after processing
        """
        outf = self._copywith(props=True, dimensions=False)
        dimlens = OrderedDict()
        for dk, dv in self.dimensions.items():
            dimlens[dk] = len(dv)

        for dk, df in dimfuncs.items():
            dv = self.dimensions[dk]
            if dk in dimfuncs:
                if verbose > 1:
                    print(dk, flush=True)
                if dk in self.variables:
                    dvar = self.variables[dk]
                    if dvar.ndim != 1:
                        dvar = np.arange(len(dv))
                else:
                    dvar = np.arange(len(dv))
                if isinstance(df, str):
                    newdl = getattr(dvar[...], df)(keepdims=True).size
                else:
                    newdl = df(dvar[:]).size
            else:
                newdl = len(dv)
            dimlens[dk] = newdl

        for dk, dv in self.dimensions.items():
            newdl = dimlens[dk]
            outf.copyDimension(dv, key=dk, dimlen=newdl)

        if verbose == 1:
            print('|' + '=' * len(self.variables.keys()) + '|', flush=True)
            print('|', end='', flush=True)

        for vark, varo in self.variables.items():
            vdims = varo.dimensions
            newvals = varo[...]
            dik = list(enumerate(vdims))
            for di, dk in dik[::-1]:
                if dk in dimfuncs:
                    if verbose > 1:
                        print(' ' * 100, '\r', vark, dk, end='\r', flush=True)
                    elif verbose > 0:
                        print('.', end='', flush=True)
                    opts = dict(axis=di, arr=newvals)
                    dfunc = dimfuncs[dk]
                    if isinstance(dfunc, dict):
                        opts.update(dfunc)
                        noopts = False
                    else:
                        opts['func1d'] = dfunc
                        noopts = True
                    if noopts and isinstance(dfunc, str):
                        newvals = getattr(newvals, dfunc)(
                            axis=di, keepdims=True)
                    else:
                        newvals = np.apply_along_axis(dfunc, di, newvals)
            newvaro = outf.copyVariable(varo, key=vark, withdata=False)
            newvaro[...] = newvals
        if verbose > 0:
            print()

        return outf

    def getTimes(self, datetype='datetime', bounds=False):
        """
        Get an array of datetime objects

        Parameters
        ----------
        datetype : string or numpy.dtype
            'datetime' or datetime64 dtype
        bounds : boolean
            get time boundaries

        Returns
        -------
        out : array
            datetime objects or array of numpy's datetype type

        Notes
        -----
        self must have a time or TFLAG variable
        """
        from PseudoNetCDF.coordutil import _parse_ref_date
        from datetime import date, datetime, timedelta, timezone
        utc = timezone.utc

        _calendaryearlike = {'noleap': 1970, '365_day': 1970,
                             'all_leap': 1972, '366_day': 1972}

        if 'time' in self.variables.keys():
            time = self.variables['time']
            timeunits = time.units.strip()
            calendar = getattr(time, 'calendar', 'gregorian').lower()
            if bounds:
                if 'time_bounds' in self.variables.keys():
                    time = self.variables['time_bounds']
                    time = np.append(time[:, 0], time[-1, 1])
                else:
                    dts = np.diff(time)
                    dt = dts.mean()
                    if (dt != dts).all():
                        warn('Time bounds are approximate')
                    time = np.append(time - dt / 2, time[-1] + dt / 2)
            if 'since' in timeunits:
                unit, base = timeunits.split(' since ')
                # Get the reference date
                refdate = _parse_ref_date(base)

                if calendar in _calendaryearlike:
                    refyear = refdate.year
                    # Get a year for relative day calculations
                    yearlike = _calendaryearlike[calendar]
                    # In that year, how many seconds and days are there
                    yearseconds = (date(yearlike + 1, 1, 1) -
                                   date(yearlike, 1, 1)).total_seconds()
                    yeardays = yearseconds / 3600 / 24

                    # Get a new reference date in yearlike
                    crefdate = datetime(yearlike, 1, 1, tzinfo=utc)
                    if refdate.month != 1 or refdate.day != 1:
                        # Get start date in yearlike
                        refcdate = datetime(
                            yearlike, refdate.month, refdate.day, tzinfo=utc)
                        # Calculate delta in years
                        addyears = (
                            crefdate - refcdate).total_seconds() / yearseconds
                    else:
                        addyears = 0
                    # Convert time to fractional years, including change in
                    # reference
                    incrdenom = {'years': 1, 'days': yeardays,
                                 'hours': yeardays * 24,
                                 'minutes': yeardays * 24 * 60,
                                 'seconds': yeardays * 24 * 60}[unit]
                    fracyearincrs = time[:] / incrdenom + addyears
                    # Split into years and days
                    yearincrs = np.array(fracyearincrs // 1).astype('i')
                    dayincrs = (fracyearincrs % 1) * yeardays
                    # Add days to the calendar year reference
                    cdays = [crefdate + timedelta(days=dayinc)
                             for dayinc in dayincrs]
                    try:
                        # Combine calendar specific month and day with new year
                        out = np.array([
                            datetime(refyear + yearinc, cday.month,
                                     cday.day, tzinfo=utc)
                            for yearinc, cday in zip(yearincrs, cdays)])
                    except Exception:
                        warn(('Years calculated from %d day year, but ' +
                              'month/days calculated for actual year. ' +
                              'Usually means data has Feb 29th in a non ' +
                              'leap year') % yeardays)
                        out = np.array([
                            datetime(refyear + yearinc, 1, 1, tzinfo=utc) +
                            timedelta(days=float(dayinc))
                            for yearinc, dayinc in zip(yearincrs, dayincrs)])

                else:
                    out = refdate + \
                        np.array([timedelta(**{unit: float(i)})
                                  for i in time[:]])
            else:
                return time
        elif 'TFLAG' in self.variables.keys():
            dates = self.variables['TFLAG'][:][:, 0, 0]
            if (dates == -635).any():
                warn('Dates of -635 set to 1970001')
                dates[dates == -635] = 1970001
            times = self.variables['TFLAG'][:][:, 0, 1]
            yyyys = (dates // 1000).astype('i')
            jjj = dates % 1000
            hours = times // 10000
            minutes = times % 10000 // 100
            seconds = times % 100
            days = jjj + (hours + minutes / 60. + seconds / 3600.) / 24.
            out = np.array([datetime(yyyy, 1, 1, tzinfo=utc) +
                            timedelta(days=day - 1)
                            for yyyy, day in zip(yyyys, days)])
            if bounds:
                if hasattr(self, 'TSTEP'):
                    tstep = getattr(self, 'TSTEP')
                    sh = tstep // 10000 * 3600
                    sm = tstep % 10000 // 100 * 60
                    ss = tstep % 100
                    dt = timedelta(seconds=sh + sm + ss)
                else:
                    dts = np.diff(out)
                    dt = dts.mean()
                    if (dt != dts).all():
                        warn('Time bounds are approximate')
                out = np.append(out, out[-1] + dt)
        elif hasattr(self, 'SDATE') and hasattr(self, 'STIME') and \
                hasattr(self, 'TSTEP') and 'TSTEP' in self.dimensions:
            jdate = self.SDATE
            if jdate < 1:
                warn(f'SDATE was {jdate}; using 1970001')
                jdate = 1970001
            hhmmss = self.STIME
            refdate = datetime.strptime(
                '%07d %06d+0000' % (jdate, hhmmss), '%Y%j %H%M%S%z')
            tstepstr = '%06d' % self.TSTEP
            timeincr = timedelta(seconds=int(tstepstr[-2:])) + \
                timedelta(minutes=int(tstepstr[-4:-2])) + \
                timedelta(hours=int(tstepstr[:-4]))
            ntimes = len(self.dimensions['TSTEP'])
            if bounds:
                ntimes += 1
            timeincrs = timeincr * np.arange(ntimes)
            out = refdate + timeincrs
        elif 'tau0' in self.variables.keys():
            out = datetime(1985, 1, 1, 0, tzinfo=utc) + \
                np.array([timedelta(hours=i)
                          for i in self.variables['tau0'][:]])
            if bounds and 'tau1' in self.variables.keys():
                oute = (
                    datetime(1985, 1, 1, 0, tzinfo=utc) +
                    np.array([
                        timedelta(hours=i)
                        for i in self.variables['tau1'][:]
                    ])
                )
                out = np.append(out, oute[-1])
        else:
            raise ValueError('cannot understand time for file')
        if datetype == 'datetime':
            return out
        else:
            return np.array(out, dtype=datetype)

    def stack(self, other, stackdim):
        """
        Concatenates all variables on stackdim

        Parameters
        ----------
        other : instance or list of PseudoNetCDFFiles
            files to add to this file along stackdim

        stackdim : str
            dimension name

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with stacked variables and dimension equal to new length
        """
        from collections.abc import Iterable
        outf = self._copywith(props=True, dimensions=False)
        if isinstance(other, Iterable):
            fs = [self] + list(other)
        else:
            fs = [self, other]
        dimensions = [f_.dimensions for f_ in fs]
        shareddims = {}
        for dimk, dim in self.dimensions.items():
            if dimk == stackdim:
                continue
            dimlens = [len(dims[dimk]) for dims in dimensions]
            if all([len(dim) == i for i in dimlens]):
                shareddims[dimk] = len(dim)
        differentdims = [set(dims.keys()).difference(
            shareddims.keys()) for dims in dimensions]
        assert(all([different.union([stackdim]) == set([stackdim])
                    for different in differentdims]))
        for dimkey in shareddims:
            ind = self.dimensions[dimkey]
            outf.copyDimension(ind, key=dimkey)
        newdl = sum([len(dims[stackdim]) for dims in dimensions])
        outf.copyDimension(self.dimensions[stackdim],
                           key=stackdim, dimlen=newdl)
        for tmpf in fs:
            for varkey, var in tmpf.variables.items():
                if stackdim not in var.dimensions:
                    if varkey in outf.variables:
                        if (
                            np.array_equal(outf.variables[varkey][...],
                                           var[...])
                        ):
                            pass
                        elif varkey not in self.dimensions:
                            warn(('Got duplicate variables for %s without ' +
                                  'stackable dimension; first value retained')
                                 % varkey)
                        continue
                    else:
                        outvals = var[...]
                else:
                    if varkey not in outf.variables.keys():
                        axisi = list(var.dimensions).index(stackdim)
                        outvals = np.ma.concatenate(
                            [f_.variables[varkey][:] for f_ in fs], axis=axisi)
                    else:
                        continue
                outvar = outf.copyVariable(var, key=varkey, withdata=False)
                outvar[...] = outvals

        return outf

    def subsetVariables(
        self, varkeys, inplace=False, exclude=False, keepcoords=True
    ):
        """
        Return a PseudoNetCDFFile with only varkeys

        Parameters
        ----------
        varkeys : iterable of strings
            keys to keep
        inplace : boolean
            if true (default false), then remove other variable from
            this file
        exclude : boolean
            if True (default False), then remove just these variables
        keepcoords : boolean
            if True (default True), keep coordinate variables

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with variables
        """
        if exclude:
            varkeys = list(set(list(self.variables)).difference(varkeys))

        varkeys = varkeys + [k for k in self.getCoords() if k not in varkeys]

        if inplace:
            outf = self
            for varkey in list(outf.variables):
                if varkey not in varkeys:
                    del outf.variables[varkey]
        else:
            outf = self._copywith(props=True, dimensions=True)
            for varkey in varkeys:
                varo = self.variables[varkey]
                outf.copyVariable(varo, key=varkey, withdata=True)
        return outf

    def sliceDimensions(self, newdims=('POINTS',), verbose=0, **dimslices):
        """
        Return a netcdflike object with dimensions sliced

        Parameters
        ----------
        newdims : iterable of strings
            names for new dimensions. When more than one iterable applies to a
            variable slice, fancy indexing removes both dimensions and creates
            a new one of the iterable lengths
        **dimslices : dictionary
            key value pairs where the key is a dimension and the value is a
            valid slice object (slices, ints or iterables) if iterables are
            provided, all iterables must be the same size and shape. If the
            arrays are not 1D, newdims must have ndim names

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with variables and dimensions sliced
        """
        outf = self._copywith(props=True, dimensions=False)
        newdimlens = OrderedDict()
        for dk, dv in self.dimensions.items():
            newdimlens[dk] = len(dv)

        isarray = {dk: not np.isscalar(dv) and not isinstance(
            dv, slice) for dk, dv in dimslices.items()}
        anyisarray = np.sum(list(isarray.values())) > 1

        if anyisarray:
            arraylens = np.array(
                [np.asarray(da).size
                 for dk, da in dimslices.items() if isarray[dk]])
            arrayshapes = np.array(
                [np.asarray(da).shape
                 for dk, da in dimslices.items() if isarray[dk]])
            arraylen = arraylens[0]
            arrayshape = arrayshapes[0]
            if (
                not (arraylens == arraylen).all() or
                not (arrayshapes == arrayshape).all()
            ):
                raise ValueError(
                    'If slicing with arrays, they must all be the same size ' +
                    'and shape')
            for dk, ia in isarray.items():
                if ia:
                    dimslices[dk] = np.asarray(dimslices[dk])

        for dk, ds in dimslices.items():
            # if anyisarray and isarray[dk]: continue
            dv = self.dimensions[dk]
            if dk in dimslices:
                if dk in self.variables:
                    dvar = self.variables[dk]
                else:
                    dvar = np.arange(len(dv))
                newdl = dvar[dimslices[dk]].size
            else:
                newdl = len(dv)
            newdimlens[dk] = newdl

        for dk, dl in newdimlens.items():
            dv = self.dimensions[dk]
            outf.copyDimension(dv, key=dk, dimlen=dl)

        if anyisarray:
            for ni, newdim in enumerate(newdims):
                outf.createDimension(newdim, arrayshape[ni])

        if verbose == 1:
            print('|' + '=' * len(self.variables.keys()) + '|', flush=True)
            print('|', end='', flush=True)
        for vark, varo in self.variables.items():
            if verbose > 1:
                print(' ' * 100, '\r', vark, end='\r', flush=True)
            elif verbose > 0:
                print('.', end='', flush=True)
            odims = vdims = varo.dimensions
            sliceo = tuple(dimslices.get(dk, slice(None)) for dk in vdims)
            isdarray = [isarray.get(dk, False) for dk in vdims]
            needsfancy = sum(isdarray) > 1
            if anyisarray and needsfancy:
                concatax = np.argmax(isdarray)
                odims = [dk for dk in vdims if not isarray.get(dk, False)]
                for newdim in newdims[::-1]:
                    odims.insert(concatax, newdim)

            newvaro = outf.copyVariable(
                varo, key=vark, dimensions=odims, withdata=False)
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
                    point_arrays.append(np.expand_dims(
                        varo[sliceoi], axis=concatax))
                newvals = np.concatenate(point_arrays, axis=concatax)
            else:
                newvals = varo[sliceo]
            try:
                newvaro[...] = newvals
            except Exception:
                newvaro[...] = newvals.reshape(newvaro.shape)

        if verbose > 0:
            print()
        return outf

    def removeSingleton(self, dimkey=None):
        """
        Return a netcdflike object with dimensions sliced

        Parameters
        ----------
        dimkey : string
            key of dimension to be evaluated for removal; if None, evaluate
            all. only singleton dimensions will be removed.

        Returns
        -------
        outf : PseudoNetCDFFile
            instance with dimensions removed
        """
        outf = self._copywith(props=True, dimensions=False)
        removed_dims = []
        for dk, d in self.dimensions.items():
            ni = len(d)
            if (dimkey is None or dk == dimkey) and ni == 1:
                removed_dims.append(dk)
            else:
                outf.copyDimension(d, key=dk, dimlen=ni)

        for vk, v in self.variables.items():
            olddims = v.dimensions
            newdims = tuple(
                [dk for dk in v.dimensions if dk not in removed_dims])
            sdims = tuple([(di, dk) for di, dk in enumerate(
                olddims) if dk not in newdims])[::-1]
            propd = dict([(pk, getattr(v, pk)) for pk in v.ncattrs()])
            ov = outf.createVariable(vk, v.dtype.char, newdims, **propd)
            outvals = v[...]
            for di, dk in sdims:
                outvals = outvals.take(0, axis=di)

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
        pncdump(self, header=True, outfile=out)
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
        new.dimensions = PseudoNetCDFDimensions()
        new._ncattrs = ()
        new._operator_exclude_vars = ()
        return new

    def __init__(self, *args, **properties):
        mode = properties.pop('mode', 'w')
        self._mode = mode
        for k, v in properties.items():
            setattr(self, k, v)

    @classmethod
    def open_mfdataset(cls, *paths, stackdim=None, **kwds):
        files = [cls(p, **kwds) for p in paths]
        file1 = files[0]
        if stackdim is None:
            for dk, dv in file1.dimensions.items():
                if dv.isunlimited():
                    stackdim = dk
            else:
                for dk in list(file1.dimensions):
                    if dk in ('TSTEP', 'time', 'Time', 't'):
                        stackdim = dk
                else:
                    raise ValueError(
                        'No dimension is unlimited or time; ' +
                        'must specify stackdim'
                    )
        return file1.stack(files[1:], stackdim=stackdim)

    def iswritable(self):
        return (self._mode[:1] in ('a', 'w') or self._mode[:2] in ('r+',))

    def __setattr__(self, k, v):
        if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
            if k not in self._ncattrs:
                self._ncattrs += (k, )
        object.__setattr__(self, k, v)

    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def setCoords(self, keys, missing='ignore'):
        """
        Set a variable as a coordinate variable

        Parameters
        ----------
        keys : iterable of strings
            keys for coord variables
        missing : string
            action if missing 'ignore', 'skip' or 'error'
            ignore - add in case used later
            skip   - do not add
            error  - raise an error

        Returns
        -------
        None

        Notes
        -----
        Coordinate variables are excluded from math
        """

        if missing == 'ignore':
            pass
        elif missing == 'skip':
            keys = [key for key in keys if key in self.variables]
        elif missing == 'error':
            invalid_keys = [key for key in keys if key not in self.variables]
            if len(invalid_keys) > 0:
                raise KeyError('File does not contain the variables: ' +
                               '{}'.format(invalid_keys))
        else:
            raise ValueError(
                'Must be ignore, skip or error; received {}'.format(missing))

        new = tuple(set(self._operator_exclude_vars).union(set(keys)))
        try:
            self._operator_exclude_vars = new
        except Exception:
            return new

    def getCoords(self):
        """
        Return a list of coordkeys
        """
        return tuple([k for k in self._operator_exclude_vars])

    def createDimension(self, name, length):
        """
        Create a dimension

        Parameters
        ----------
        name : string
            name for dimension
        length : integer
            maximum length of dimension

        Returns
        -------
        dim : PseudoNetCDFDimensions
            new dimension
        """
        dim = self.dimensions[name] = PseudoNetCDFDimension(self, name, length)
        return dim

    def copyDimension(self, dim, key=None, dimlen=None, unlimited=None):
        if dimlen is None:
            dimlen = len(dim)

        if key is None:
            key = dim.name

        if unlimited is None:
            unlimited = dim.isunlimited()

        if isinstance(self, netcdf):
            if unlimited:
                ndv = self.createDimension(key, None)
            else:
                ndv = self.createDimension(key, dimlen)
        else:
            ndv = self.createDimension(key, dimlen)
            ndv.setunlimited(unlimited)

        return ndv

    def copyVariable(self, var, key=None, dtype=None, dimensions=None,
                     fill_value=None, withdata=True):
        """
        Copy var into self as vark

        Parameters
        ----------
        var : PseudoNetCDFVariable
            netCDF4.Variable-like object (must have ncattrs and setncatts)
        key : string
            key for variable in self (can be omitted if var has name,
            standard_name, or long_name)
        dtype : string or numpy.dtype
            change the data type to dtype
        dimensions : iterable of strings
            change the dimensions to dimensions
        fill_value : integer or flaot
            change the fill_value to this values
        withdata : boolean
            default True, copies data

        Returns
        -------
        myvar : PseuodNetCDFVairable
            copy of var in this file
        """
        if key is None:
            for propk in ['_name', 'name', 'standard_name', 'long_name']:
                if hasattr(var, propk):
                    key = getattr(var, propk)
                    break
            else:
                raise AttributeError(
                    'varkey must be supplied because var has no name, ' +
                    'standard_name or long_name')

        if withdata:
            try:
                vals = var[:]
            except Exception:
                vals = var[...]
        else:
            vals = var

        if dtype is None:
            dtype = vals.dtype

        if dimensions is None:
            dimensions = var.dimensions
        if fill_value is None:
            for pk in ('fill_value', 'missing_value', '_FillValue'):
                fill_value = getattr(var, pk, None)
                if fill_value is not None:
                    break

        myvar = self.createVariable(
            key, dtype, dimensions, fill_value=fill_value)
        attrs = OrderedDict()
        for propk in var.ncattrs():
            attrs[propk] = getattr(var, propk)
        myvar.setncatts(attrs)
        if withdata:
            try:
                myvar[:] = vals[:]
            except Exception:
                myvar[...] = vals[...]
        return myvar

    def createVariable(self, name, type, dimensions, fill_value=None,
                       **properties):
        """
        Create a variable

        Parameters
        ----------
        name : string
            name for new variable
        type : string or numpy dtype
            code (e.g., 'f', 'i', 'd')
        dimensions : tuple of strigns
            dimension keys that can be found in objects' dimensions dictionary

        Returns
        -------
        var : new variable
        """
        import numpy as np
        varopt = self.get_varopt()
        if fill_value is None:
            fill_value = varopt.get('fill_value', None)

        for pk, pv in varopt.items():
            if pk != 'fill_value':
                properties.setdefault(pk, pv)

        if fill_value is None:
            for pk in 'missing_value _FillValue'.split():
                fill_value = properties.get(pk, None)
                if fill_value is not None:
                    break
        else:
            properties['fill_value'] = fill_value

        if type == 'S':
            type = 'c'
        if (
            isinstance(properties.get('values', 1), np.ma.MaskedArray) or
            fill_value is not None
        ):
            var = self.variables[name] = PseudoNetCDFMaskedVariable(
                self, name, type, dimensions, **properties)
        else:
            var = self.variables[name] = PseudoNetCDFVariable(
                self, name, type, dimensions, **properties)
        return var

    def close(self):
        """
        Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
        """
        pass

    def dump(self, *args, **kwds):
        from PseudoNetCDF.pncdump import pncdump
        pncdump(self, *args, **kwds)

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
        return pncbo(op='+', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __sub__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='-', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __mul__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='*', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __truediv__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='/', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __floordiv__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='//', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __pow__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='**', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __and__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='&', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __or__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='|', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __xor__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='^', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __mod__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='%', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __lt__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='<', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __gt__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='>', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __eq__(self, lhs):
        if isinstance(lhs, (NetCDFFile, PseudoNetCDFFile)):
            from ._functions import pncbo
            return pncbo(op=' == ', ifile1=self, ifile2=lhs, verbose=0,
                         coordkeys=self._operator_exclude_vars)
        else:
            return lhs.__eq__(self)

    def __le__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='<=', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __ge__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='>=', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    def __ne__(self, lhs):
        from ._functions import pncbo
        return pncbo(op='!=', ifile1=self, ifile2=lhs, verbose=0,
                     coordkeys=self._operator_exclude_vars)

    slice = sliceDimensions
    apply = applyAlongDimensions
    subset = subsetVariables
    sync = close
    flush = close


_ncvarkwds = [
    'varname', 'datatype', 'dimensions', 'zlib', 'complevel', 'shuffle',
    'fletcher32', 'contiguous', 'chunksizes', 'endian',
    'least_significant_digit', 'fill_value'
]


class netcdf(PseudoNetCDFFile, NetCDFFile):
    def createDimension(self, *args, **kwds):
        return NetCDFFile.createDimension(self, *args, **kwds)

    def createVariable(self, *args, **kwds):
        for pk, pv in self.get_varopt().items():
            kwds.setdefault(pk, pv)
        propkwds = OrderedDict()
        for pk in set(kwds).difference(_ncvarkwds):
            propkwds[pk] = kwds.pop(pk)

        outval = NetCDFFile.createVariable(self, *args, **kwds)
        for pk, pv in propkwds.items():
            outval.setncattr(pk, pv)

        return outval

    def setncattr(self, k, v):
        return NetCDFFile.setncattr(self, k, v)

    def __setattr__(self, k, v):
        return NetCDFFile.__setattr__(self, k, v)

    def __delattr__(self, k):
        NetCDFFile.__delattr__(self, k)

    def __new__(cls, *args, **kwds):
        return NetCDFFile.__new__(cls, *args, **kwds)

    def setCoords(self, keys, missing='ignore'):
        new = PseudoNetCDFFile.setCoords(self, keys, missing=missing)
        self.__dict__['_operator_exclude_vars'] = tuple(set(new))

    def __init__(self, *args, **kwds):
        NetCDFFile.__init__(self, *args, **kwds)
        self.__dict__['_mode'] = kwds.get('mode', 'r')
        self.__dict__['_operator_exclude_vars'] = ()
        self.__dict__['_varopt'] = {}
        coords = [k for k in self.dimensions if k in self.variables]
        self.setCoords(coords)

    @property
    def _mode(self):
        return self.__dict__.get('_mode', 'r')

    def get_dest(self):
        """
        Returns the path where a new file is created on some action

        If None, a file is created in memory.
        Else, a netcdf file is created on disk.
        """
        return self.__dict__.get('_destination', None)

    def set_dest(self, path, **options):
        """
        Sets the path where a new file is created on some action

        Arguments
        ---------
        path : path for new file
        options : options for new file creation

        Returns
        -------
        None
        """
        options['filename'] = path
        self.__dict__['_destination'] = options

    def get_varopt(self):
        """
        Get options

        Arguments
        ---------
        None

        Returns
        -------
        options : dictionary of optiosn
        """
        return self.__dict__.get('_varopt', {}).copy()

    def set_varopt(self, **options):
        """
        Set options to be used when creating any Variable

        Arguments
        ---------
        options : options for new Variable creation

        Returns
        -------
        None
        """
        self.__dict__['_varopt'] = options

    def ncattrs(self):
        return NetCDFFile.ncattrs(self)

    @classmethod
    def from_ncf(cls, infile):
        outf = PseudoNetCDFFile()
        for pk in infile.ncattrs():
            pv = getattr(infile, pk)
            setattr(outf, pk, pv)

        for dk, dv in infile.dimensions.items():
            outf.copyDimension(dv, key=dk)

        for vk, vv in infile.variables.items():
            outf.copyVariable(vv, key=vk)

        return outf

    def _newlike(self):
        """
        Internal function to return a file of the same class if a
        PsueoNetCDFFile
        """
        if self.get_dest() is not None:
            outf = netcdf(**self.get_dest())
        else:
            outf = PseudoNetCDFFile()
        outf.set_varopt(**self.get_varopt())
        return outf

    @classmethod
    def isMine(cls, path, *args, **kwds):
        """
        True if this file or object can be identified
        for use by this class. Useful to override for
        classes that can be initialized from disk.
        """
        tmpf = open(path, 'rb')
        tmpf.seek(0, 0)
        cdftest = tmpf.read(3)
        tmpf.seek(1, 0)
        hdftest = tmpf.read(3)
        tmpf.close()
        if cdftest == b'CDF':
            return True
        elif hdftest == b'HDF':
            try:
                cls(path, *args, **kwds)
                return True
            except Exception:
                return False
        else:
            return False

    def close(self):
        try:
            return NetCDFFile.close(self)
        except Exception as e:
            warn(str(e))

    def __del__(self):
        self.close()


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
        func: function
            must take a key and provides a PseudoNetCDFVariable
        keys: iterable of strings
            keys that the dictionary should act as if it has
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
        if k not in self.__keys:
            self.__keys.append(k)

    def keys(self):
        return tuple(
            self.__keys +
            [k for k in dict.keys(self) if k not in self.__keys])

    def __iter__(self):
        for k in self.keys():
            yield k

    def __len__(self):
        return len([k for k in self.keys()])

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def __contains__(self, k):
        return k in [k for k in self.keys()]


if __name__ == '__main__':
    unittest.main()
