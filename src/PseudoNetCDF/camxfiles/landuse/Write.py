import numpy as np


def ncf2landuse(ncffile, outpath):
    nland = len(ncffile.dimensions['LANDUSE'])
    nrows = len(ncffile.dimensions['ROW'])
    ncols = len(ncffile.dimensions['COL'])
    newstyle = getattr(ncffile, '_newstyle', True)
    if newstyle:
        _fland_dtype = np.dtype(dict(names = ['SPAD1', 'KEY', 'EPAD1', 'SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '8>S', '>i', '>i', '(%d, %d, %d)>f' % (nland, nrows, ncols), '>i']))
        _other_dtype = np.dtype(dict(names = ['SPAD1', 'KEY', 'EPAD1', 'SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '8>S', '>i', '>i', '(%d, %d)>f' % (nrows, ncols), '>i']))
    else:
        _fland_dtype = np.dtype(dict(names = ['SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '(%d, %d, %d)>f' % (nland, nrows, ncols), '>i']))
        _other_dtype = np.dtype(dict(names = ['SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '(%d, %d)>f' % (nrows, ncols), '>i']))

    outfile=open(outpath,'wb')
    keys = [key for key in ['FLAND', 'VAR1', 'LAI', 'TOPO'] if key in
 ncffile.variables.keys()]
    keyandvar = [(key, np.empty(shape = (1,), dtype = {'FLAND': _fland_dtype}.get(key, _other_dtype))) for key in keys]
    for key, var in keyandvar:
        invar = ncffile.variables[key]
        if newstyle:
            if key == 'FLAND'.strip():
                key = 'LUCAT%02d' % nland
        
            var['SPAD1'] = 8
            var['KEY'] = key.ljust(8)
            var['EPAD1'] = 8

        var['SPAD2'] = invar.size * 4
        var['DATA'][:] = invar[:]
        var['EPAD2'] = invar.size * 4
        var.tofile(outfile)

    outfile.flush()
    return outfile

