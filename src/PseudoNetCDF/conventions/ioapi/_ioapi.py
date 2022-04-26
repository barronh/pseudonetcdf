from PseudoNetCDF.pncwarn import warn
import numpy as np
_withlatlon = False
try:
    import pyproj
    _withlatlon = True
except Exception:
    warn('pyproj could not be found, so IO/API coordinates cannot be ' +
         'converted to lat/lon; to fix, install pyproj or basemap ' +
         '(e.g., `pip install pyproj)`')


def add_lay_coordinates(ifileo):
    if 'LAY' in ifileo.dimensions:
        nlay = len(ifileo.dimensions['LAY'])
    elif hasattr(ifileo, 'NLAYS'):
        nlay = ifileo.NLAYS
    elif hasattr(ifileo, 'VGLVLS'):
        nlay = len(ifileo.VGLVLS) - 1
    else:
        return
    if 'layer' not in ifileo.variables.keys():
        var = ifileo.createVariable('layer', 'd', ('LAY',))
        var[:] = np.arange(nlay, dtype='d')
        var.units = 'model layers'
        var.standard_name = 'layer'
    if 'level' not in ifileo.variables.keys():
        var = ifileo.createVariable('level', 'd', ('LAY',))
        if hasattr(ifileo, 'VGLVLS'):
            var[:] = (ifileo.VGLVLS[:-1] + ifileo.VGLVLS[1:]) / 2
        else:
            var[:] = np.arange(nlay, dtype='d')
        var.units = 'sigma'
        var.positive = 'down'
        var.standard_name = 'level'
        var._CoordinateAxisType = "GeoZ"
        var._CoordinateZisPositive = "down"


def add_time_variables(ifileo):
    add_time_variable(ifileo, 'time')
    add_time_variable(ifileo, 'time_bounds')


def add_time_variable(ifileo, key):
    from datetime import datetime, timedelta, timezone
    strptime = datetime.strptime

    sdate = int(abs(ifileo.SDATE))
    if sdate < 1400000:
        sdate += 2000000
    sdate = datetime(sdate // 1000, 1, 1, tzinfo=timezone.utc) + \
        timedelta(days=(sdate % 1000) - 1)
    sdate = sdate + timedelta(days=ifileo.STIME / 240000.)
    rdate = strptime('1970-01-01 00:00:00+0000', '%Y-%m-%d %H:%M:%S%z')
    if ifileo.TSTEP == 0:
        tmpseconds = 0
    else:
        tmp = ('%06d' % ifileo.TSTEP)
        htmp = tmp[:2]
        mtmp = tmp[2:4]
        stmp = tmp[4:]
        tmpseconds = 3600 * int(htmp) + 60 * int(mtmp) + int(stmp)

    time_unit = "seconds since %s" % (rdate.strftime('%Y-%m-%d %H:%M:%S%z'),)
    if 'TFLAG' in ifileo.variables:
        tflag = ifileo.variables['TFLAG'][:, 0]
        jdays = tflag[:, 0]
        hhmmsss = tflag[:, 1]
        if jdays[0] == 0:
            jdays = [ifileo.SDATE]
        dates = np.array([strptime('%7d %06d+0000' % (jday, hhmmss),
                                   '%Y%j %H%M%S%z')
                          for jday, hhmmss in zip(jdays, hhmmsss)])
        time = np.array([(date - rdate).total_seconds() for date in dates])
        if key == 'time_bounds':
            time = np.array([time, time + tmpseconds]).T
            dims = ('TSTEP', 'tnv')
            if 'tnv' not in ifileo.dimensions.keys():
                ifileo.createDimension('tnv', 2)
        else:
            dims = ('TSTEP',)
    else:
        sdate = strptime('%07d %06d+0000' %
                         (getattr(ifileo, 'SDATE', 1970001),
                          getattr(ifileo, 'STIME', 0)),
                         '%Y%j %H%M%S%z')
        off = (sdate - rdate).total_seconds()
        if key == 'time':
            tmax = max(1, len(ifileo.dimensions['TSTEP']))
            time = np.arange(0, tmax, dtype='i') * tmpseconds + off
            dims = ('TSTEP',)
        elif key == 'time_bounds':
            tmax = max(1, len(ifileo.dimensions['TSTEP'])) + 1
            time = np.arange(0, tmax, dtype='i') * tmpseconds + off
            time = time.repeat(2, 0)[1:-1].reshape(-1, 2)
            dims = ('TSTEP', 'tnv')
            if 'tnv' not in ifileo.dimensions.keys():
                ifileo.createDimension('tnv', 2)
        else:
            raise KeyError(
                'time variables are time and time_bounds, got %s' % key)
    if key not in ifileo.variables.keys():
        var = ifileo.createVariable(key, time.dtype.char, dims)
    else:
        var = ifileo.variables[key]
    var[:] = time

    var.units = time_unit
    if key == 'time':
        var._CoordinateAxisType = "Time"
        var.bounds = 'time_bounds'

    var.long_name = ("synthesized time coordinate from SDATE, " +
                     "STIME, STEP global attributes")


# 0 Clarke 1866
# 1 Clarke 1880
# 2 Bessel
# 3 New International 1967
# 4 International 1909
# 5 WGS 72
# 6 Everest
# 7 WGS 66
# 8 GRS 1980
# 9 Airy
# 10 Modified Everest
# 11 Modified Airy
# 12 WGS 84
# 13 Southeast Asia
# 14 Australian National
# 15 Krassovsky
# 16 Hough
# 17 Mercury 1960
# 18 Modified Mercury 1968
# 19 Normal Sphere
# 20 MM5 Sphere
# 21 WRF-NMM Sphere
_AXIS = np.array([6378206.4, 6378249.145, 6377397.155, 6378157.5, 6378388.0,
                  6378135.0, 6377276.3452, 6378145.0, 6378137.0, 6377563.396,
                  6377304.063, 6377340.189, 6378137.0, 6378155., 6378160.0,
                  6378245.0, 6378270.0, 6378166.0, 6378150.0, 6370997.0,
                  6370000.0, 6371200.0])
_BXIS = np.array([6356583.8, 6356514.86955, 6356078.96284, 6356772.2,
                  6356911.94613, 6356750.519915, 6356075.4133, 6356759.769356,
                  6356752.314140, 6356256.91, 6356103.039, 6356034.448,
                  6356752.314245, 6356773.3205, 6356774.719, 6356863.0188,
                  6356794.343479, 6356784.283666, 6356768.337303, 6370997.0,
                  6370000.0, 6371200.0])


def get_ioapi_sphere():
    import os
    ENV_IOAPI_ISPH = os.environ.get('IOAPI_ISPH', None)
    if ENV_IOAPI_ISPH is None:
        ENV_IOAPI_ISPH = '6370000.'
        warn('IOAPI_ISPH is assumed to be ' +
             ENV_IOAPI_ISPH + '; consistent with WRF')
    isph_parts = [eval(ip) for ip in ENV_IOAPI_ISPH.split(' ')]
    if len(isph_parts) > 2:
        raise ValueError('IOAPI_ISPH must be 1 or 2 parameters (got: %s)' %
                         str(isph_parts))
    elif len(isph_parts) == 2:
        return isph_parts
    elif isph_parts[0] >= 0 and isph_parts[0] < _AXIS.size:
        return _AXIS[isph_parts[0]], _BXIS[isph_parts[0]]
    else:
        return isph_parts * 2


_gdnames = {1: "latitude_longitude", 2: "lambert_conformal_conic",
            7: "mercator", 6: "polar_stereographic"}


def getmapdef(ifileo, add=True):
    gridname = _gdnames[ifileo.GDTYP]
    if add:
        mapdef = ifileo.createVariable(gridname, 'i', ())
    else:
        from PseudoNetCDF import PseudoNetCDFVariable
        mapdef = PseudoNetCDFVariable(ifileo, gridname, 'i', ())
    mapdef.grid_mapping_name = gridname
    if mapdef.grid_mapping_name == "latitude_longitude":
        mapdef.latitude_of_projection_origin = ifileo.YORIG
        mapdef.longitude_of_central_meridian = ifileo.XORIG
    elif mapdef.grid_mapping_name == "polar_stereographic":
        mapdef.latitude_of_projection_origin = ifileo.P_ALP * 90
        mapdef.straight_vertical_longitude_from_pole = ifileo.P_GAM
        mapdef.standard_parallel = np.array([ifileo.P_BET], dtype='d')
        # mapdef.standard_parallel = np.array([ifileo.P_BET], dtype = 'f')
    else:
        mapdef.standard_parallel = np.array(
            [ifileo.P_ALP, ifileo.P_BET], dtype='d')
    if mapdef.grid_mapping_name != "latitude_longitude":
        mapdef.latitude_of_projection_origin = ifileo.YCENT
        mapdef.longitude_of_central_meridian = ifileo.XCENT
        mapdef.false_northing = -ifileo.YORIG
        mapdef.false_easting = -ifileo.XORIG
        mapdef._CoordinateTransformType = "Projection"
        mapdef._CoordinateAxes = "x y"
    ioapi_sphere = get_ioapi_sphere()
    mapdef.semi_major_axis = getattr(ifileo, 'semi_major_axis', getattr(
        ifileo, 'earth_radius', ioapi_sphere[0]))
    mapdef.semi_minor_axis = getattr(ifileo, 'semi_minor_axis', getattr(
        ifileo, 'earth_radius', ioapi_sphere[1]))
    return mapdef


def add_lcc_coordinates(ifileo):
    mapdef = getmapdef(ifileo)
    gridname = mapdef.grid_mapping_name
    if 'PERIM' in ifileo.dimensions.keys():
        xdim = 'PERIM'
        ydim = 'PERIM'
        _x = np.arange(-ifileo.XCELL, (ifileo.NCOLS + 1) * ifileo.XCELL,
                       ifileo.XCELL) + ifileo.XORIG + ifileo.XCELL / 2.
        _y = np.arange(-ifileo.YCELL, (ifileo.NROWS + 1) * ifileo.YCELL,
                       ifileo.YCELL) + ifileo.YORIG + ifileo.YCELL / 2.
        bx = _x[1:]
        by = _y[0].repeat(ifileo.NCOLS + 1)
        ex = _x[-1].repeat(ifileo.NROWS + 1)
        ey = _y[1:]
        tx = _x[0:-1]
        ty = _y[-1].repeat(ifileo.NCOLS + 1)
        wx = _x[0].repeat(ifileo.NROWS + 1)
        wy = _y[:-1]
        lcc_x = x = np.concatenate([bx, ex, tx, wx])
        lcc_y = y = np.concatenate([by, ey, ty, wy])
        XCELL = ifileo.XCELL
        YCELL = ifileo.YCELL
        lcc_xe = np.array([x - XCELL / 2., x + XCELL /
                           2., x + XCELL / 2., x - XCELL / 2.]).T
        lcc_ye = np.array([y - YCELL / 2., y - YCELL /
                           2., y + YCELL / 2., y + YCELL / 2.]).T
        latlon_dim = ('PERIM',)
        latlone_dim = ('PERIM', 'nv')
        latlon_coord = 'PERIM'
    else:
        xdim = 'COL'
        ydim = 'ROW'
        latlon_dim = (ydim, xdim)
        latlon_coord = 'latitude longitude'
        latlone_dim = (ydim, xdim, 'nv')
        xe = np.arange(0, ifileo.NCOLS) * ifileo.XCELL + ifileo.XORIG
        ye = np.arange(0, ifileo.NROWS) * ifileo.YCELL + ifileo.YORIG
        lcc_xe, lcc_ye = np.meshgrid(xe, ye)
        lcc_xe = np.concatenate([lcc_xe[:, :, None],
                                 lcc_xe[:, :, None] + ifileo.XCELL,
                                 lcc_xe[:, :, None] + ifileo.XCELL,
                                 lcc_xe[:, :, None]], axis=2)
        lcc_ye = np.concatenate([lcc_ye[:, :, None], lcc_ye[:, :, None],
                                 lcc_ye[:, :, None] + ifileo.YCELL,
                                 lcc_ye[:, :, None] + ifileo.YCELL], axis=2)
        x = np.arange(0, ifileo.NCOLS) * ifileo.XCELL + \
            ifileo.XCELL / 2. + ifileo.XORIG
        y = np.arange(0, ifileo.NROWS) * ifileo.YCELL + \
            ifileo.YCELL / 2. + ifileo.YORIG
        lcc_x, lcc_y = np.meshgrid(x, y)

    if _withlatlon:
        if ifileo.GDTYP == 2:
            mapstrs = ['+proj=lcc',
                       '+a=%s' % mapdef.semi_major_axis,
                       '+b=%s' % mapdef.semi_minor_axis,
                       '+lon_0=%s' % mapdef.longitude_of_central_meridian,
                       '+lat_1=%s' % mapdef.standard_parallel[0],
                       '+lat_2=%s' % mapdef.standard_parallel[1],
                       '+lat_0=%s' % mapdef.latitude_of_projection_origin]
            mapstr = ' '.join(mapstrs)
            mapproj = pyproj.Proj(mapstr)
        elif ifileo.GDTYP == 6:
            mapstr = ('+proj=stere +a={3} +b={4} ' +
                      '+lon_0={0} +lat_0={1} +lat_ts={2}').format(
                          mapdef.straight_vertical_longitude_from_pole,
                          mapdef.latitude_of_projection_origin,
                          mapdef.standard_parallel[0],
                          mapdef.semi_major_axis,
                          mapdef.semi_minor_axis)
            mapproj = pyproj.Proj(mapstr)
        elif ifileo.GDTYP == 7:
            mapstr = '+proj=merc +a=%s +b=%s +lat_ts=0 +lon_0=%s' % (
                mapdef.semi_major_axis, mapdef.semi_minor_axis,
                mapdef.longitude_of_central_meridian)
            mapproj = pyproj.Proj(mapstr)
        elif ifileo.GDTYP == 1:
            def mapproj(x, y, inverse):
                return (x, y)
        lon, lat = mapproj(lcc_x, lcc_y, inverse=True)
        lone, late = mapproj(lcc_xe.ravel(), lcc_ye.ravel(), inverse=True)
        lone = lone.reshape(*lcc_xe.shape)
        late = late.reshape(*lcc_ye.shape)

    if 'x' not in ifileo.variables.keys() and xdim in ifileo.dimensions:
        """
        Not necessary for cdo
        """
        var = ifileo.createVariable('x', x.dtype.char, (xdim,))
        var[:] = x[:]
        var.units = 'km'
        var._CoordinateAxisType = "GeoX"
        var.long_name = ("synthesized coordinate from XORIG XCELL " +
                         "global attributes")

    if 'y' not in ifileo.variables.keys() and ydim in ifileo.dimensions:
        """
        Not necessary for cdo
        """
        var = ifileo.createVariable('y', x.dtype.char, (ydim,))
        var[:] = y[:]
        var.units = 'km'
        var._CoordinateAxisType = "GeoY"
        var.long_name = ("synthesized coordinate from YORIG YCELL " +
                         "global attributes")

    if _withlatlon and 'latitude' not in ifileo.variables.keys():
        var = ifileo.createVariable('latitude', lat.dtype.char, latlon_dim)
        var[:] = lat
        var.units = 'degrees_north'
        var.standard_name = 'latitude'
        var.bounds = 'latitude_bounds'
        var.coordinates = latlon_coord

    if _withlatlon and 'longitude' not in ifileo.variables.keys():
        var = ifileo.createVariable('longitude', lon.dtype.char, latlon_dim)
        var[:] = lon
        var.units = 'degrees_east'
        var.standard_name = 'longitude'
        var.bounds = 'longitude_bounds'
        var.coordinates = latlon_coord

    if _withlatlon:
        for dk, dl in zip(latlone_dim, late.shape):
            if dk not in ifileo.dimensions:
                ifileo.createDimension(dk, dl)

    if _withlatlon and 'latitude_bounds' not in ifileo.variables.keys():
        var = ifileo.createVariable(
            'latitude_bounds', lat.dtype.char, latlone_dim)
        var[:] = late
        var.units = 'degrees_north'
        var.standard_name = 'latitude_bounds'

    if _withlatlon and 'longitude_bounds' not in ifileo.variables.keys():
        var = ifileo.createVariable(
            'longitude_bounds', lon.dtype.char, latlone_dim)
        var[:] = lone
        var.units = 'degrees_east'
        var.standard_name = 'longitude_bounds'

    for varkey in ifileo.variables.keys():
        var = ifileo.variables[varkey]
        # this must have been a fix for dictionaries that
        # reproduced variables on demand
        # we should find a better fix for this
        # try:
        #    ifileo.variables[varkey] = var
        # except Exception:
        #    pass
        olddims = list(var.dimensions)
        dims = [d for d in olddims]  # Why was I excluding time  if d != 'time'
        if _withlatlon:
            def io2cf(x):
                return {'ROW': 'latitude', 'COL': 'longitude',
                        'TSTEP': 'time', 'LAY': 'level'}.get(x, x)
            dims = [io2cf(d) for d in olddims]

        if olddims != dims:
            if (varkey not in ('latitude', 'longitude') and
                ('PERIM' in dims or
                 ('latitude' in dims and 'longitude' in dims))):
                try:
                    var.coordinates = ' '.join(dims)
                    var.grid_mapping = gridname
                except Exception as e:
                    warn(('coordinates="{0}" and gridmapping="{1}" ' +
                          'not added to variables:\n\t{3}'
                          ).format(' '.join(dims), gridname, varkey, e),
                         category=UserWarning)


def add_ioapi_from_cf(ifileo, coordkeys=[]):
    from ...coordutil import gettimes
    times = gettimes(ifileo)
    jdays = np.array([int(t.strftime('%Y%j')) for t in times])
    itimes = np.array([int(t.strftime('%H%M%S')) for t in times])
    outkeys = [k.ljust(16) for k in ifileo.variables.keys()
               if k not in ('ETFLAG', 'TFLAG') and k not in coordkeys]
    setattr(ifileo, 'VAR-LIST', ''.join(outkeys))
    setattr(ifileo, 'NVARS', len(outkeys))
    ifileo.createDimension('VAR', ifileo.NVARS)
    setattr(ifileo, 'SDATE', jdays[0])
    setattr(ifileo, 'STIME', itimes[0])
    setattr(ifileo, 'EDATE', jdays[-1])
    setattr(ifileo, 'ETIME', itimes[-1])
    tflag = ifileo.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag[:, :, 0] = jdays[:, None]
    tflag[:, :, 1] = itimes[:, None]
    tflag.units = "<YYYYDDD,HHMMSS>"
    tflag.long_name = "TFLAG           "
    tflag.var_desc = ("Timestep-valid flags: " +
                      " (1) YYYYDDD or (2) HHMMSS").ljust(80)


add_ioapi_from_ioapi = add_ioapi_from_cf


def add_cf_from_ioapi(ifileo, coordkeys=[]):
    try:
        add_lay_coordinates(ifileo)
    except Exception:
        pass
    try:
        add_time_variables(ifileo)
    except Exception:
        pass
    if ifileo.GDTYP in _gdnames:
        add_lcc_coordinates(ifileo)
    else:
        raise TypeError('IOAPI is only aware of LCC (GDTYP=2)')
    ifileo.Conventions = 'CF-1.6'
