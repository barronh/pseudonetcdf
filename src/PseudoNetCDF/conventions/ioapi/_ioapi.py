from warnings import warn
import numpy as np
_withlatlon = False
for impstmt in ['import pyproj', 'from mpl_toolkits.basemap import pyproj']:
    try:
        exec(impstmt)
    except ImportError:
        pass
    else:
        _withlatlon = True
        break
else:
    warn('pyproj could not be found, so IO/API coordinates cannot be converted to lat/lon; to fix, install pyproj or basemap (e.g., `pip install pyproj)`')

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
        var[:] = np.arange(nlay, dtype = 'd')
        var.units = 'model layers'
        var.standard_name = 'layer';
    if 'level' not in ifileo.variables.keys():
        var = ifileo.createVariable('level', 'd', ('LAY',))
        var[:] = np.arange(nlay, dtype = 'd')
        var.units = 'sigma'
        var.positive = 'down'
        var.standard_name = 'level';
        var._CoordinateAxisType = "GeoZ" ;
        var._CoordinateZisPositive = "down";

def add_time_variables(ifileo):
    add_time_variable(ifileo, 'time')
    add_time_variable(ifileo, 'time_bounds')

def add_time_variable(ifileo, key):
    from datetime import datetime, timedelta
    sdate = int(abs(ifileo.SDATE))
    if sdate < 1400000:
        sdate += 2000000
    sdate = datetime(sdate // 1000, 1, 1) + timedelta(days = (sdate % 1000) - 1)
    sdate = sdate + timedelta(days = ifileo.STIME / 240000.)
    if ifileo.TSTEP == 0:
        tmpseconds = 0;
    else:
        tmp = ('%06d' % ifileo.TSTEP)
        htmp = tmp[:2]
        mtmp = tmp[2:4]
        stmp = tmp[4:]
        tmpseconds = 3600 * int(htmp) + 60 * int(mtmp) + int(stmp)

    time_unit = "seconds since %s UTC" % (sdate.strftime('%Y-%m-%d %H:%M:%S'),)
    if 'TFLAG' in ifileo.variables:
        tflag = ifileo.variables['TFLAG'][:, 0]
        jdays = tflag[:, 0]
        hhmmsss = tflag[:, 1]
        if jdays[0] == 0:
            jdays = [ifileo.SDATE]
        dates = np.array([datetime.strptime('%7d %06d' % (jday, hhmmss), '%Y%j %H%M%S') for jday, hhmmss in zip(jdays, hhmmsss)])
        time = np.array([(date - sdate).total_seconds() for date in dates])
        if key == 'time_bounds':
            time = np.array([time, time + tmpseconds]).T
        
    else:
        if key == 'time':
            time = np.arange(0, max(1, len(ifileo.dimensions['TSTEP'])), dtype = 'i') * tmpseconds
            dims = ('TSTEP',)
        elif key == 'time_bounds':
            time = np.arange(0, max(1, len(ifileo.dimensions['TSTEP'])) + 1, dtype = 'i') * tmpseconds
            time = time.repeat(2, 0)[1:-1].reshape(-1, 2)
            dims = ('TSTEP','tnv')
            if 'tnv' not in ifileo.dimensions.keys():
                ifileo.createDimension('tnv', 2)
        else:
            raise KeyError('time variables are time and time_bounds, got %s' % key)
    if key not in ifileo.variables.keys():
        var = ifileo.createVariable(key, time.dtype.char, dims)
    else:
        var = ifileo.variables[key]
    var[:] = time
    var.units = time_unit
    if key == 'time':
        var._CoordinateAxisType = "Time" ;
        var.bounds = 'time_bounds'
    
    var.long_name = "synthesized time coordinate from SDATE, STIME, STEP global attributes" ;

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
_AXIS = np.array([6378206.4,6378249.145,6377397.155,6378157.5,6378388.0,6378135.0,6377276.3452,6378145.0,6378137.0,6377563.396,6377304.063,6377340.189,6378137.0,6378155.,6378160.0,6378245.0,6378270.0,6378166.0,6378150.0,6370997.0,6370000.0,6371200.0])
_BXIS = np.array([6356583.8,6356514.86955,6356078.96284,6356772.2,6356911.94613,6356750.519915,6356075.4133,6356759.769356,6356752.314140,6356256.91,6356103.039,6356034.448,6356752.314245,6356773.3205,6356774.719,6356863.0188,6356794.343479,6356784.283666,6356768.337303,6370997.0,6370000.0,6371200.0])
def get_ioapi_sphere():
    import os
    isph_parts = map(eval, os.environ.get('IOAPI_ISPH', '6370000.').split(' '))
    if len(isph_parts) > 2:
        raise ValueError('IOAPI_ISPH must be 1 or 2 parameters (got: %s)' % str(isph_parts))
    elif len(isph_parts) == 2:
        return isph_parts
    elif isph_parts > 0 and isph_parts < _AXIS.size :
        return _AXIS[isph_parts], _BXIS[isph_parts]
    else:
        return isph_parts * 2

def add_lcc_coordinates(ifileo, lccname = 'LambertConformalProjection'):
    mapdef = ifileo.createVariable(lccname, 'i', ())
    if ifileo.GDTYP == 2:
        mapdef.grid_mapping_name = "lambert_conformal_conic" ;
    elif ifileo.GDTYP == 7:
        mapdef.grid_mapping_name = "equatorial_mercator" ;
    mapdef.standard_parallel = np.array([ifileo.P_ALP, ifileo.P_BET], dtype = 'f')
    ioapi_sphere = get_ioapi_sphere()
    mapdef.semi_major_axis = getattr(ifileo, 'semi_major_axis', getattr(ifileo, 'earth_radius', ioapi_sphere[0]))
    mapdef.semi_minor_axis = getattr(ifileo, 'semi_minor_axis', getattr(ifileo, 'earth_radius', ioapi_sphere[1]))
    mapdef.latitude_of_projection_origin = ifileo.YCENT
    mapdef.longitude_of_central_meridian = ifileo.XCENT
    mapdef._CoordinateTransformType = "Projection" ;
    mapdef._CoordinateAxes = "x y" ;

    if 'PERIM' in ifileo.dimensions.keys():
        xdim = 'PERIM'
        ydim = 'PERIM'
        _x = np.arange(-ifileo.XCELL, (ifileo.NCOLS + 1) * ifileo.XCELL, ifileo.XCELL) + ifileo.XORIG + ifileo.XCELL / 2.
        _y = np.arange(-ifileo.YCELL, (ifileo.NROWS + 1) * ifileo.YCELL, ifileo.YCELL) + ifileo.YORIG + ifileo.YCELL / 2.
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
        lcc_xe = np.array([x - ifileo.XCELL / 2., x + ifileo.XCELL / 2., x + ifileo.XCELL / 2., x - ifileo.XCELL / 2.]).T
        lcc_ye = np.array([y - ifileo.YCELL / 2., y - ifileo.YCELL / 2., y + ifileo.YCELL / 2., y + ifileo.YCELL / 2.]).T
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
        lcc_xe = np.concatenate([lcc_xe[:, :, None], lcc_xe[:, :, None] + ifileo.XCELL, lcc_xe[:, :, None] + ifileo.XCELL, lcc_xe[:, :, None]], axis = 2)
        lcc_ye = np.concatenate([lcc_ye[:, :, None], lcc_ye[:, :, None], lcc_ye[:, :, None] + ifileo.YCELL, lcc_ye[:, :, None] + ifileo.YCELL], axis = 2)
        x = np.arange(0, ifileo.NCOLS) * ifileo.XCELL + ifileo.XCELL / 2. + ifileo.XORIG
        y = np.arange(0, ifileo.NROWS) * ifileo.YCELL  + ifileo.YCELL / 2. + ifileo.YORIG
        lcc_x, lcc_y = np.meshgrid(x, y)

    if _withlatlon:
        if ifileo.GDTYP == 2:
            mapstr = '+proj=lcc +lon_0=%s +lat_1=%s +lat_2=%s +a=%s +b=%s +lat_0=%s' % (mapdef.longitude_of_central_meridian, mapdef.standard_parallel[0], mapdef.standard_parallel[1], mapdef.semi_major_axis, mapdef.semi_minor_axis, mapdef.latitude_of_projection_origin,) 
            mapproj = pyproj.Proj(mapstr)
        elif ifileo.GDTYP == 7:
            mapstr = '+proj=merc +a=%s +b=%s +lat_ts=0 +lon_0=%s' % (mapdef.semi_major_axis, mapdef.semi_minor_axis, mapdef.longitude_of_central_meridian)
            mapproj = pyproj.Proj(mapstr)
        elif ifileo.GDTYP == 1:
            mapproj = lambda x, y, inverse: (x, y)
        lon, lat = mapproj(lcc_x, lcc_y, inverse = True)
        lone, late = mapproj(lcc_xe, lcc_ye, inverse = True)
        
    if 'x' not in ifileo.variables.keys():
        """
        Not necessary for cdo
        """
        var = ifileo.createVariable('x', x.dtype.char, (xdim,))
        var[:] = x[:]
        var.units = 'km'
        var._CoordinateAxisType = "GeoX" ;
        var.long_name = "synthesized coordinate from XORIG XCELL global attributes" ;

    
    if 'y' not in ifileo.variables.keys():
        """
        Not necessary for cdo
        """
        var = ifileo.createVariable('y', x.dtype.char, (ydim,))
        var[:] = y[:]
        var.units = 'km'
        var._CoordinateAxisType = "GeoY" ;
        var.long_name = "synthesized coordinate from YORIG YCELL global attributes" ;


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
        var.standard_name = 'longitude';
        var.bounds = 'longitude_bounds'
        var.coordinates = latlon_coord

    if _withlatlon:
        for dk, dl in zip(latlone_dim, late.shape):
            if not dk in ifileo.dimensions:
                ifileo.createDimension(dk, dl)

    if _withlatlon and 'latitude_bounds' not in ifileo.variables.keys():
        var = ifileo.createVariable('latitude_bounds', lat.dtype.char, latlone_dim)
        var[:] = late
        var.units = 'degrees_north'
        var.standard_name = 'latitude_bounds'

    if _withlatlon and 'longitude_bounds' not in ifileo.variables.keys():
        var = ifileo.createVariable('longitude_bounds', lon.dtype.char, latlone_dim)
        var[:] = lone
        var.units = 'degrees_east'
        var.standard_name = 'longitude_bounds';

    for varkey in ifileo.variables.keys():
        var = ifileo.variables[varkey]
        try:
            ifileo.variables[varkey] = var
        except:
            pass
        olddims = list(var.dimensions)
        if _withlatlon:
            dims = map(lambda x: {'ROW': 'latitude', 'COL': 'longitude', 'TSTEP': 'time', 'LAY': 'level'}.get(x, x), olddims)
        dims = [d for d in dims] # Why was I excluding time  if d != 'time'
        if olddims != dims:
            if ('PERIM' in dims or 
                ('latitude' in dims and 'longitude' in dims)
               ) and varkey not in ('latitude', 'longitude'):
                var.coordinates = ' '.join(dims)
                var.grid_mapping = lccname

def add_cf_from_ioapi(ifileo):
    try:
        add_lay_coordinates(ifileo)
    except:
        pass
    try:
        add_time_variables(ifileo)
    except:
        pass
    if ifileo.GDTYP in (1, 2, 7):
        add_lcc_coordinates(ifileo)
    else:
        raise TypeError('IOAPI is only aware of LCC (GDTYP=2)')
    ifileo.Conventions = 'CF-1.6'