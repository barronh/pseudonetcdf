from collections import OrderedDict
import numpy as np
from .. import PseudoNetCDFFile, PseudoNetCDFVariables
from datetime import datetime, timedelta
import re


_shortnamere = re.compile('(.+)-(.+)')


def _parseforecastseconds(md):
    fs = md['GRIB_FORECAST_SECONDS']
    fs = fs.replace(' sec', '').replace(' s', '')
    return float(fs)


def _parseshortname(md):
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
    else:
        try:
            level = float(level)
        except Exception:
            coordname = shortname
            level = 0

    return coordname, level


class grib2(PseudoNetCDFFile):
    def __init__(self, path):
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
        for i in range(1, nr + 1):
            mesg = ds.GetRasterBand(i)
            md = mesg.GetMetadata()
            coordname, level = _parseshortname(md)
            zcoords.setdefault(coordname, set())
            zcoords[coordname].add(level)
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
        ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)
        for zk, zvals in zcoords.items():
            zv = self.createVariable(zk, 'd', (zk,))
            zv.units = zk
            zv[:] = sorted(zvals)
        xv = self.createVariable('x', 'd', ('x',))
        xv[:] = np.linspace(ulx, lrx, ds.RasterXSize)
        yv = self.createVariable('y', 'd', ('y',))
        yv[:] = np.linspace(uly, lry, ds.RasterYSize)
        key = keys[0]
        tmpvar = self.variables[key]
        for dk, dv in zip(tmpvar.dimensions, tmpvar.shape):
            self.createDimension(dk, dv)

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

        if clip:
            i = np.minimum(np.maximum(0, i), nx - 1)
            j = np.minimum(np.maximum(0, j), ny - 1)
        if mask:
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
            zk, level = _parseshortname(md)
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
