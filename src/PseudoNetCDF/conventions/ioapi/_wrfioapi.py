import numpy as np
_coorddict = dict(west_east = 'longitude', south_north = 'latitude', Time = 'time', bottom_top = 'altitude',
                  west_east_stag = 'longitude', south_north_stag = 'latitude', Time_stag = 'time', bottom_top_stag = 'altitude',)
def add_cf_from_wrfioapi(ifile):
    for invark, outvark in [('XLONG', 'longitude'), ('XLAT', 'latitude')]:
        invar = ifile.variables[invark]
        outvar = ifile.createVariable(outvark, invar.dtype.char, invar.dimensions[1:])
        for pk in invar.ncattrs():
            setattr(outvar, pk, getattr(invar, pk))
        outvar[:] = invar[0]
        
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
    ifile.Conventions = 'CF-1.6'
    