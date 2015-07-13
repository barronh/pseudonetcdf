import numpy as np
_coorddict = dict(west_east = 'longitude', south_north = 'latitude', Time = 'time', bottom_top = 'altitude',
                  west_east_stag = 'longitude', south_north_stag = 'latitude', Time_stag = 'time', bottom_top_stag = 'altitude',)
def add_cf_from_wrfioapi(ifile):
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
    except: pass
            
    for k in ifile.variables.keys():
        var = ifile.variables[k]
        var.coordinates = ' '.join([_coorddict.get(dk, dk) for dk in var.dimensions])
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
    for varkey in ifile.variables.keys():
        var = ifile.variables[varkey]
        try:
            ifile.variables[varkey] = var
        except:
            pass
        olddims = list(var.dimensions)
        dims = map(lambda x: {'west_east': 'latitude', 'north_south': 'longitude', 'Time': 'time', 'LAY': 'level'}.get(x, x), olddims)
        dims = [d for d in dims] # Why was I excluding time  if d != 'time'
        if olddims != dims:
            if ('PERIM' in dims or 
                ('latitude' in dims and 'longitude' in dims)
               ) and varkey not in ('latitude', 'longitude'):
                var.coordinates = ' '.join(dims)

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
    projname = {1: "lcc", 2: 'npstere', 3: 'merc'}
    
    if projname == 'merc':
        var = ifile.createVariable('mercator', 'i', ())
        var.grid_mapping_name = 'mercator'
        var.longitude_of_projection_origin = ifile.STAND_LON
        var.standard_parallel = ifile.TRUELAT1
        var.false_easting = len(ifile.dimensions['west_east']) / 2
        var.false_northing = len(ifile.dimensions['south_north']) / 2
    elif projname == "lcc":
        var = ifile.createVariable('mercator', 'i', ())
        var.grid_mapping_name = 'lambert_conformal_conic'
        var.standard_parallel = np.array([ifile.TRUELAT1, ifile.TRUELAT1])
        var.longitude_of_central_meridian = ifile.STAND_LON
        var.latitude_of_projection_origin = ifile.REF_LAT
        var.false_easting = 0
        var.false_northing = 0

    elif projname == 'npstere':
        if ifile.TRUELAT1 > 0:
            projname = 'npstere'
        else:
            projname = 'spstere'
        var = ifile.createVariable('mercator', 'i', ())
        var.grid_mapping_name = 'polar_stereographic'
        var.straight_vertical_longitude_from_pole = ifile.STAND_LON
        var.latitude_of_projection_origin = ifile.TRUELAT1
        var.standard_parallel = ifile.REF_LALT
        var.false_easting = 0
        var.false_northing = 0

