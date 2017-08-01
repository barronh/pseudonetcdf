import numpy as np
_coorddict = dict(west_east = 'longitude', south_north = 'latitude', Time = 'time', bottom_top = 'altitude',
                  west_east_stag = 'longitude', south_north_stag = 'latitude', Time_stag = 'time', bottom_top_stag = 'altitude',)
def add_cf_from_wrfioapi(ifile, coordkeys = []):
    try:
        for invark, outvark in [('XLONG', 'longitude'), ('XLAT', 'latitude')]:
            try:
                invar = ifile.variables[invark]
            except KeyError:
                invark += '_M'
                invar = ifile.variables[invark]
            outvar = ifile.createVariable(outvark, invar.dtype.char, invar.dimensions[1:])
            for pk in invar.ncattrs():
                setattr(outvar, pk, getattr(invar, pk))
            outvar[:] = invar[0]
    except Exception as e:
        print(e)
        pass
    
    # Add time coordinate if not available.
    try:
        outvar = ifile.createVariable('time', 'i', ('Time',))
        invar = ifile.variables['Times']
        for pk in invar.ncattrs():
            setattr(outvar, pk, getattr(invar, pk))
        outvar.units = 'seconds since 1985-01-01T00:00:00Z'
        invals = invar[:].copy().view('S19')
        sdate = np.datetime64('1985-01-01T00:00:00Z')
        timesince = (np.array([np.datetime64(time[0].replace('_', 'T') + 'Z') for time in invals]) - sdate).astype('l')
        outvar[:] = timesince
    except: pass
    # add projection
    grid_mapping_name = get_proj(ifile)
    # Add x and y coordinates
    try:
        gridvar = ifile.variables[grid_mapping_name]
        xvar = ifile.createVariable('x', 'd', ('west_east',))
        xvar.units = 'm'
        xvar.standard_name = 'x'
        xvar[:] = np.arange(xvar.size, dtype = 'd')*ifile.DX-gridvar.false_easting  + ifile.DX/2
        yvar = ifile.createVariable('y', 'd', ('south_north',))
        yvar.units = 'm'
        yvar.standard_name = 'y'
        yvar[:] = np.arange(yvar.size, dtype = 'd')*ifile.DY-gridvar.false_northing  + ifile.DY/2
    except Exception as e:
        raise e
    for k in ifile.variables.keys():
        var = ifile.variables[k]
        olddims = [dk for dk in var.dimensions]
        newdims = [_coorddict.get(dk, dk) for dk in var.dimensions]
        if olddims != newdims:
            var.coordinates = ' '.join(newdims)
            var.grid_mapping = grid_mapping_name
    

    ifile.Conventions = 'CF-1.6'
    
def get_proj(ifile):
    """
    MAP_PROJ - Model projection [1=Lambert, 2=polar stereographic, 3=mercator, 6=lat-lon]  (required)
    TRUELAT1 - required for MAP_PROJ = 1, 2, 3 (defaults to 0 otherwise)
    TRUELAT2 - required for MAP_PROJ = 6 (defaults to 0 otherwise)
    STAND_LON - Standard longitude used in model projection (required)
    REF_LON, REF_LON - A reference longitude and latitude (required)
    KNOWNI, KNOWNJ - The I and J locations of REF_LON and REF_LAT (required)
    POLE_LAT - optional for MAP_PROJ = 6 (defaults to 90 otherwise)
    POLE_LAT - optional for MAP_PROJ = 6 (defaults to 0 otherwise)
    DX, DY - required for MAP_PROJ = 1, 2, 3 (defaults to 0 otherwise)
    LATINC, LONINC - required for MAP_PROJ = 6 (defaults to 0 otherwise)
    """
    
    projname = {1: "lambert_conformal_conic", 2: 'polar_stereographic', 3: 'mercator'}[ifile.MAP_PROJ]
    var = ifile.createVariable(projname, 'i', ())
    var.grid_mapping_name = projname
    var.earth_radius = 6370000.
    var.false_easting = 0
    var.false_northing = 0
    x0_from_cen = len(ifile.dimensions['west_east']) / 2 * ifile.DX
    y0_from_cen = len(ifile.dimensions['south_north']) / 2 * ifile.DY
    if projname == 'mercator':
        var.longitude_of_projection_origin = ifile.STAND_LON
        var.standard_parallel = ifile.TRUELAT1
    elif projname == "lambert_conformal_conic":
        var.standard_parallel = np.array([ifile.TRUELAT1, ifile.TRUELAT2])
        var.longitude_of_central_meridian = ifile.STAND_LON
        var.latitude_of_projection_origin = ifile.MOAD_CEN_LAT
    elif projname == 'polar_stereographic':
        var.straight_vertical_longitude_from_pole = ifile.STAND_LON
        var.latitude_of_projection_origin = ifile.TRUELAT1
        var.standard_parallel = ifile.MOAD_CEN_LAT
    
    from PseudoNetCDF.coordutil import getproj4_from_cf_var
    from mpl_toolkits.basemap import pyproj
    projstr = getproj4_from_cf_var(var)
    proj = pyproj.Proj(projstr)
    gcx, gcy = proj(ifile.CEN_LON, ifile.CEN_LAT)
    glx = gcx - x0_from_cen
    gly = gcy - y0_from_cen
    var.false_easting = -glx
    var.false_northing = -gly
    return projname
