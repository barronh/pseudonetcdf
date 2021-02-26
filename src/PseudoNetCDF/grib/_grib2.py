from collections import OrderedDict
import numpy as np
from .. import PseudoNetCDFFile, PseudoNetCDFVariables
from datetime import datetime, timedelta
import re


_shortnamere = re.compile('(.+)-(.+)')


def _parseforecastseconds(md):
    """
    Convenience function to parse time in seconds from a raster
    GetMetadata_Dict

    Arguments
    ---------
    md : dict
        Dictionary representing the raster metadata.

    Returns
    -------
    fs : float
        Forecast time in seconds
    """
    fs = md['GRIB_FORECAST_SECONDS']
    fs = fs.replace(' sec', '').replace(' s', '')
    return float(fs)


def _parseshortname(md):
    """
    Convenience function to parse coordinate name, level, and bounds from
    a raster GetMetadata_Dict result.

    Arguments
    ---------
    md : dict
        Dictionary representing the raster metadata.

    Returns
    -------
    coordname, level, bounds
    """
    shortname = md['GRIB_SHORT_NAME']
    m = _shortnamere.match(shortname)
    if m is None:
        coordname = shortname
        level = '0'
    else:
        level, coordname = m.groups()

    if '-' in level[1:]:
        low, hi = level.split('-')
        mid = (float(low) + float(hi)) / 2
        level = mid
        bounds = tuple(sorted([low, hi]))
        coordname = coordname + f'_{low}-{hi}'
    else:
        try:
            level = float(level)
            bounds = None
        except Exception:
            coordname = shortname
            level = 0
            bounds = None

    return coordname, level, bounds


class grib2(PseudoNetCDFFile):
    def __init__(self, path):
        """
        Arguments
        ---------
        path : str
            path to a grib2 file

        Returns
        -------
        None

        Notes
        -----
        Initializes PseudoNetCDFFile that uses raster bands to construct 4D
        data.
          * Times are based on parsing the GRIB_FORECAST_SECONDS, and
          * Levels are based on parsing the GRIB_SHORT_NAME
            * Short names are parsed to split the z-coordinate system (e.g.,
              HTGL, ISBL, SFC, etc) from the numeric descriptor.
            * Numeric descriptors are used to construct coordinates
          * Variables are the combination of GRIB_ELEMENT and the z-coordinate
            system. this allows for multiple variables (e.g., TMP_SFC and
            TMP_HTGL). Each variable will have all bands with that coordinate
            type.
        """
        try:
            from osgeo import gdal
            from osgeo import osr
        except ImportError:
            raise IOError('Reading grib2 files requires gdal library')

        PseudoNetCDFFile.__init__(self)
        ds = self._dataset = gdal.Open(path, gdal.GA_ReadOnly)
        self.wktstring = ds.GetProjection()
        srs = osr.SpatialReference()
        srs.ImportFromWkt(self.wktstring)
        self.proj4string = srs.ExportToProj4()
        nr = self._nrast = ds.RasterCount
        k2r = self._key2raster = OrderedDict()
        k2u = self._key2unit = OrderedDict()
        times = set()
        zcoords = OrderedDict()
        zbounds = OrderedDict()
        for i in range(1, nr + 1):
            mesg = ds.GetRasterBand(i)
            md = mesg.GetMetadata()
            coordname, level, lbounds = _parseshortname(md)
            zcoords.setdefault(coordname, set())
            zcoords[coordname].add(level)
            if lbounds is not None:
                zbounds.setdefault(coordname, set())
                zbounds[coordname].add(lbounds)
            key = md['GRIB_ELEMENT'] + '_' + coordname
            unit = md['GRIB_UNIT']
            k2r.setdefault(key, [])
            k2u.setdefault(key, [])
            k2r[key].append(i)
            k2u[key].append(unit)
            times.add(_parseforecastseconds(md))

        self.createDimension('time', len(times))
        for zk, zvals in zcoords.items():
            self.createDimension(zk, len(zvals))

        self.createDimension('x', ds.RasterXSize)
        self.createDimension('y', ds.RasterYSize)
        self.createDimension('nv', 2)
        self.reftime = md['GRIB_REF_TIME'].strip()
        assert(self.reftime.endswith('sec UTC'))
        refseconds = int(self.reftime.split()[0])
        refdate = datetime(1970, 1, 1) + timedelta(seconds=refseconds)
        keys = list(k2r)
        self.variables = PseudoNetCDFVariables(
            self._getvar, keys + ['time', 'x', 'y']
        )
        tv = self.createVariable('time', 'd', ('time',))
        tv.long_name = 'time'
        tv.description = 'GRIB_FORECAST_SECONDS'
        tv.units = refdate.strftime('seconds since %Y-%m-%dT%H:%M:%S+0000')
        tv[:] = sorted(times)
        for zk, zvals in zcoords.items():
            zv = self.createVariable(zk, 'd', (zk,))
            zv.units = zk
            zv[:] = sorted(zvals)

        for zk, zbnds in zbounds.items():
            try:
                zv = self.createVariable(zk + '_bounds', 'd', (zk, 'nv'))
                zv.units = zk
                zv[:] = list(sorted(zbnds))
            except Exception:
                import warnings
                warnings.warn(f'Unable to add {zk} bounds values {zbnds}')

        xv = self.createVariable('x', 'd', ('x',))
        xv.bounds = 'x_bounds'
        self.createVariable('x_bounds', 'd', ('x', 'nv'))
        yv = self.createVariable('y', 'd', ('y',))
        yv.bounds = 'y_bounds'
        self.createVariable('y_bounds', 'd', ('y', 'nv'))
        ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
        self.reset_proj(
            ulx=ulx, xres=xres, xskew=xskew,
            uly=uly, yres=yres, yskey=yskew
        )

    def plot(self, *args, plottype='x-y', **kwds):
        """Thin wrapper around PseudoNetCDFFile to set default plottype to x-y
        For docs see PseudoNetCDFFile.plot
        """
        return PseudoNetCDFFile.plot(self, *args, **kwds)

    def reset_proj(self, proj4string=None, wrf=False, **kwds):
        """
        Overwrite the coordinate variables x/y to change the projected
        coordinates.

        Arguments
        ---------
        proj4string : str
            A new projection definition. It may change any part.
        wrf : bool
            If True, then recalculate the urx and ury using a symmetric WRF
            domain based on lat_0 and lon_0 as the origin of counting.
        kwds : keywords
            Potential key words include ulx, uly, xres, yres, lrx, or lry
            Other keywords that are passed will be ignored.

        Returns
        -------
        None

        Notes
        -----
        This will rarely be used.

        It was added due to previous experience where the projection did not
        accurately represent the data, which was on a WRF domain. The grib
        radius was approximate, but could be improved to match WRF.
        """
        from osgeo import osr
        kwds = kwds.copy()

        xv = self.variables['x']
        xbv = self.variables['x_bounds']
        yv = self.variables['y']
        ybv = self.variables['y_bounds']
        xres = (xbv[:, 1] - xbv[:, 0]).mean(0)
        yres = (ybv[:, 1] - ybv[:, 0]).mean(0)
        ulx = xbv[0, 0]
        uly = ybv[0, 0]
        kwds.setdefault('xres', xres)
        kwds.setdefault('yres', yres)
        kwds.setdefault('ulx', ulx)
        kwds.setdefault('uly', uly)

        nx = xv.size
        ny = yv.size

        if proj4string is not None:
            self.proj4string = proj4string
            srs = osr.SpatialReference()
            srs.ImportFromProj4(self.proj4string)
            self.wkstring = srs.ExportToWkt()
        if wrf:
            x0 = -(nx + 1) / 2
            y0 = -(ny + 1) / 2
            ulx, uly = x0 * kwds['xres'], y0 * kwds['yres']
            print(ulx, uly)
            kwds.setdefault('ulx', ulx)
            kwds.setdefault('uly', uly)

        lrx = kwds['ulx'] + (nx * kwds['xres'])
        lry = kwds['uly'] + (ny * kwds['yres'])
        kwds.setdefault('lrx', lrx)
        kwds.setdefault('lry', lry)
        xedges = np.linspace(kwds['ulx'], kwds['lrx'], nx + 1)
        xbv[:, 0] = xedges[:-1]
        xbv[:, 1] = xedges[1:]
        xv[:] = xbv.mean(1)
        yedges = np.linspace(kwds['uly'], kwds['lry'], ny + 1)
        ybv[:, 0] = yedges[:-1]
        ybv[:, 1] = yedges[1:]
        yv[:] = ybv.mean(1)

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
        proj = self.getproj()
        x, y = proj(lon, lat)
        i = self.val2idx('x', x)
        yp = self.variables['y'][:]
        if np.diff(yp).mean() < 0:
            fp = np.arange(yp.size)[::-1] + 0.5
            jf = np.interp(y, yp[::-1], fp, left=yp.size, right=0)
            j = jf.round(0).astype('i')
        else:
            j = self.val2idx('y', y)

        nx = len(self.dimensions['x'])
        ny = len(self.dimensions['y'])

        if clean == 'clip':
            i = np.minimum(np.maximum(0, i), nx - 1)
            j = np.minimum(np.maximum(0, j), ny - 1)
        if clean == 'mask':
            i = np.ma.masked_greater(np.ma.masked_less(i, 0), nx - 1)
            j = np.ma.masked_greater(np.ma.masked_less(j, 0), ny - 1)

        return i, j

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
        x = self.variables['x'][i]
        y = self.variables['y'][j]
        return self.xy2ll(x, y)

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
        if withgrid:
            raise ValueError('grib2 getproj with grid')

        if projformat == 'proj4':
            return self.proj4string
        elif projformat != 'pyproj':
            raise ValueError(
                'grib2 getproj only supports projformat="pyproj"'
                + f'; got {projformat}'
            )

        import pyproj
        pstr = self.proj4string
        proj = pyproj.Proj(pstr, preserve_units=withgrid)

        return proj

    def _getvar(self, key):
        """
        Internal function used to create variables on the fly.

        Arguments
        ---------
        key : str
            Grib Raster band key, usually GRIB_ELEMENT hyphen and the last
            part of GRIB_SHORT_NAME that cooresponds to the z-coordinate name

        Returns
        -------
        pvar : PseudoNetCDFVariable
            A PseudoNetCDFVariable with metadata and values read from the
            RasterBands
        """
        k2u = self._key2unit
        myunits = k2u[key]
        assert(all([myunits[0] == u for u in myunits]))
        zdim = '_'.join(key.split('_')[1:])
        out = self.createVariable(
            key, 'f', ('time', zdim, 'y', 'x'), fill_value=-1e20
        )
        unit = myunits[0]
        if unit.startswith('[') and unit.endswith(']'):
            unit = unit[1:-1]
        out.units = unit
        out.standard_name = key
        out[:] = np.ma.masked
        times = self.variables['time'][:]
        levels = None
        for idx in self._key2raster[key]:
            mesg = self._dataset.GetRasterBand(idx)
            md = mesg.GetMetadata_Dict()
            zk, level, lbounds = _parseshortname(md)
            if levels is None:
                levels = self.variables[zk]

            time = _parseforecastseconds(md)
            data_array = mesg.ReadAsArray()
            li = np.where(levels == level)[0]
            ti = np.where(times == time)[0]
            out[ti, li] = data_array

        return out


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import pycno

    f = grib2('gdas.t00z.pgrb2.1p00.f000')
    x = f.variables['x'] - 0.5
    y = f.variables['y'] - 0.5
    p = plt.pcolormesh(x, y, f.variables['TMP-0-SFC'][0, 0])
    cno = pycno.cno()
    cno.getfeatures('MWDB_Coasts_Countries_3.cnob')
    cf = cno._cachedfeatures['MWDB_Coasts_Countries_3.cnob']
    for i in range(len(cf)):
        cf[i] = (cf[i][0] % 360, cf[i][1])

    cno.draw(ax=p.axes)
    p.axes.figure.savefig('test.png')
