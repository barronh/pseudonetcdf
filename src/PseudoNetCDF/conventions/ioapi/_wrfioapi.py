_coorddict = dict(west_east = 'longitude', south_north = 'latitude', Time = 'time', bottom_top = 'altitude',
                  west_east_stag = 'longitude', south_north_stag = 'latitude', Time_stag = 'time', bottom_top_stag = 'altitude',)
def add_cf_from_wrfioapi(ifile):
    for invark, outvark in [('XLONG', 'longitude'), ('XLAT', 'latitude')]:
        invar = ifile.variables[invark]
        outvar = ifile.createVariable(outvark, invar.dtype.char, invar.dimensions[1:])
        for pk in invar.ncattrs():
            setattr(outvar, pk, getattr(invar, pk))
    for k in ifile.variables.keys():
        var = ifile.variables[k]
        var.coordinates = ' '.join([_coorddict.get(dk, dk) for dk in var.dimensions])