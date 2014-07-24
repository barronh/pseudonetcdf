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
    warn('pyproj could not be found, so IO/API coordinates cannot be converted to lat/lon')

def add_time_variables(ifileo):
    add_time_variable(ifileo, 'time')
    add_time_variable(ifileo, 'time_bounds')

def add_time_variable(ifileo, key):
    if key not in ifileo.variables.keys():
        from datetime import datetime, timedelta
        sdate = int(abs(ifileo.SDATE))
        if sdate < 1400000:
            sdate += 2000000
        sdate = datetime(sdate // 1000, 1, 1) + timedelta(days = (sdate % 1000) - 1)
        if ifileo.TSTEP == 0:
            tmpseconds = 0;
        else:
            tmp = ('%06d' % ifileo.TSTEP)
            htmp = tmp[:2]
            mtmp = tmp[2:4]
            stmp = tmp[4:]
            tmpseconds = 3600 * int(htmp) + 60 * int(mtmp) + int(stmp)
    
        time_unit = "seconds since %s 00:00:00 UTC" % (sdate.strftime('%Y-%m-%d'),)
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
        var = ifileo.createVariable(key, time.dtype.char, dims)
        var[:] = time
        var.units = time_unit
        if key == 'time':
            var._CoordinateAxisType = "Time" ;
            var.bounds = 'time_bounds'
        
        var.long_name = "synthesized time coordinate from SDATE, STIME, STEP global attributes" ;
        

def add_lcc_coordinates(ifileo, lccname = 'LambertConformalProjection'):
    mapdef = ifileo.createVariable(lccname, 'i', ())
    if ifileo.GDTYP == 2:
        mapdef.grid_mapping_name = "lambert_conformal_conic" ;
    elif ifileo.GDTYP == 7:
        mapdef.grid_mapping_name = "equatorial_mercator" ;
    mapdef.standard_parallel = np.array([ifileo.P_ALP, ifileo.P_BET], dtype = 'f')
    mapdef.earth_radius = 6371000. # Add search for earth radius later
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
            mapstr = '+proj=lcc +lon_0=%s +lat_1=%s +lat_2=%s +a=%s +b=%s +lat_0=%s' % (mapdef.longitude_of_central_meridian, mapdef.standard_parallel[0], mapdef.standard_parallel[1], mapdef.earth_radius, mapdef.earth_radius, mapdef.latitude_of_projection_origin,) 
        elif ifileo.GDTYP == 7:
            mapstr = '+proj=merc +a=%s +b=%s +lat_ts=0 +lon_0=%s' % (mapdef.earth_radius, mapdef.earth_radius, mapdef.longitude_of_central_meridian)
        mapproj = pyproj.Proj(mapstr)
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

    if 'layer' not in ifileo.variables.keys() and 'LAY' in ifileo.dimensions:
        var = ifileo.createVariable('layer', lon.dtype.char, ('LAY',))
        var[:] = np.arange(len(ifileo.dimensions['LAY']))
        var.units = 'model layers'
        var.standard_name = 'layer';

    if 'level' not in ifileo.variables.keys() and 'LAY' in ifileo.dimensions:
        var = ifileo.createVariable('level', lon.dtype.char, ('LAY',))
        var[:] = np.arange(len(ifileo.dimensions['LAY']))
        var.units = 'sigma'
        var.positive = 'down'
        var.standard_name = 'level';
        var._CoordinateAxisType = "GeoZ" ;
        var._CoordinateZisPositive = "down";

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
    add_time_variables(ifileo)
    if ifileo.GDTYP in (2, 7):
        add_lcc_coordinates(ifileo)
    else:
        raise TypeError('IOAPI is only aware of LCC (GDTYP=2)')
    ifileo.Conventions = 'CF-1.6'