import numpy as np
try:
    from mpl_toolkits.basemap import pyproj
    _withlatlon = True
except:
    _withlatlon = False

def add_time_variable(ifileo):
    if 'time' not in ifileo.variables.keys():
        from datetime import datetime, timedelta
        sdate = int(abs(ifileo.SDATE))
        if sdate < 1400000:
            sdate += 2000000
        sdate = datetime(sdate // 1000, 1, 1) + timedelta(days = (sdate % 1000) - 1)
        if ifileo.TSTEP == 0:
            tmpseconds = 0;
        else:
            tmp = ('%06s' % ifileo.TSTEP)
            htmp = tmp[:2]
            mtmp = tmp[2:4]
            stmp = tmp[4:]
            tmpseconds = 3600 * int(htmp) + 60 * int(mtmp) + int(stmp)
    
        time_unit = "seconds since %s 00:00:00 UTC" % (sdate.strftime('%Y-%m-%d'),)
        time = np.arange(0, len(ifileo.dimensions['TSTEP'])) * tmpseconds
        var = ifileo.createVariable('time', time.dtype.char, ('TSTEP',))
        var[:] = time
        var.units = time_unit
        var._CoordinateAxisType = "Time" ;
        var.long_name = "synthesized time coordinate from SDATE, STIME, STEP global attributes" ;

def add_lcc_coordinates(ifileo, lccname = 'LambertConformalProjection'):
    lccdef = ifileo.createVariable(lccname, 'i', ())
    lccdef.grid_mapping_name = "lambert_conformal_conic" ;
    lccdef.latitude_of_projection_origin = ifileo.YCENT
    lccdef.longitude_of_central_meridian = ifileo.XCENT
    lccdef.standard_parallel = np.array([ifileo.P_ALP, ifileo.P_BET], dtype = 'f')
    lccdef.earth_radius = 6371000. # Add search for earth radius later
    lccdef._CoordinateTransformType = "Projection" ;
    lccdef._CoordinateAxes = "x y" ;

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
        latlon_dim = ('PERIM',)
        latlone_dim = ('PERIM_STAG',)
    else:
        xdim = 'COL'
        ydim = 'ROW'
        latlon_dim = (ydim, xdim)
        latlone_dim = (ydim + '_STAG', xdim + '_STAG')
        xe = np.arange(0, (ifileo.NCOLS + 1) * ifileo.XCELL, ifileo.XCELL) + ifileo.XORIG
        ye = np.arange(0, (ifileo.NROWS + 1) * ifileo.YCELL, ifileo.YCELL) + ifileo.YORIG + ifileo.YCELL
        lcc_xe, lcc_ye = np.meshgrid(xe, ye)
        x = np.arange(0, ifileo.NCOLS * ifileo.XCELL, ifileo.XCELL) + ifileo.XORIG + ifileo.XCELL / 2.
        y = np.arange(0, ifileo.NROWS * ifileo.YCELL, ifileo.YCELL) + ifileo.YORIG + ifileo.YCELL / 2.
        lcc_x, lcc_y = np.meshgrid(x, y)

    if _withlatlon: lcc = pyproj.Proj('+proj=lcc +lon_0=%s +lat_1=%s +lat_2=%s +a=%s +lat_0=%s' % (lccdef.longitude_of_central_meridian, lccdef.standard_parallel[0], lccdef.standard_parallel[1], lccdef.earth_radius, lccdef.latitude_of_projection_origin,)  )


    lon, lat = lcc(lcc_x, lcc_y, inverse = True)
    lone, late = lcc(lcc_xe, lcc_ye, inverse = True)
        
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

    if _withlatlon and 'longitude' not in ifileo.variables.keys():
        var = ifileo.createVariable('longitude', lon.dtype.char, latlon_dim)
        var[:] = lon
        var.units = 'degrees_east'
        var.standard_name = 'longitude';

    if _withlatlon and latlone_dim[0] not in ifileo.dimensions.keys():
        ifileo.createDimension(latlone_dim[0], len(ifileo.dimensions[latlon_dim[0]]) + 1)
    if _withlatlon and latlone_dim[1] not in ifileo.dimensions.keys():
        ifileo.createDimension(latlone_dim[1], len(ifileo.dimensions[latlon_dim[1]]) + 1)

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

    if 'layer' not in ifileo.variables.keys():
        var = ifileo.createVariable('layer', lon.dtype.char, ('LAY',))
        var[:] = np.arange(len(ifileo.dimensions['LAY']))
        var.units = 'model layers'
        var.standard_name = 'layer';

    if 'level' not in ifileo.variables.keys():
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
        dims = [d for d in dims if d != 'time']
        if olddims != dims:
            if ('PERIM' in dims or 
                ('latitude' in dims and 'longitude' in dims)
               ) and varkey not in ('latitude', 'longitude'):
                var.coordinates = ' '.join(dims)
                var.grid_mapping = lccname

def add_cf_from_ioapi(ifileo):
    add_time_variable(ifileo)
    if ifileo.GDTYP == 2:
        add_lcc_coordinates(ifileo)
    else:
        raise TypeError('IOAPI is only aware of LCC (GDTYP=2)')
