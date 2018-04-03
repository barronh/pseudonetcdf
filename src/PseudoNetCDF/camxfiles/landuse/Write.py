from __future__ import unicode_literals
import numpy as np
from PseudoNetCDF._getwriter import registerwriter


def ncf2landuse(ncffile, outpath):
    nland = len(ncffile.dimensions['LANDUSE'])
    nrows = len(ncffile.dimensions['ROW'])
    ncols = len(ncffile.dimensions['COL'])
    newstyle = getattr(ncffile, '_newstyle', True)
    dt3dfmt = '(%d, %d, %d)>f' % (nland, nrows, ncols)
    dt2dfmt = '(%d, %d)>f' % (nrows, ncols)
    if newstyle:
        _fland_dtype = np.dtype(dict(names=['SPAD1', 'KEY', 'EPAD1', 'SPAD2',
                                            'DATA', 'EPAD2'],
                                     formats=['>i', '8>S', '>i', '>i', dt3dfmt,
                                              '>i']))
        _other_dtype = np.dtype(dict(names=['SPAD1', 'KEY', 'EPAD1', 'SPAD2',
                                            'DATA', 'EPAD2'],
                                     formats=['>i', '8>S', '>i', '>i', dt2dfmt,
                                              '>i']))
    else:
        _fland_dtype = np.dtype(dict(names=['SPAD2', 'DATA', 'EPAD2'],
                                     formats=['>i', dt3dfmt, '>i']))
        _other_dtype = np.dtype(dict(names=['SPAD2', 'DATA', 'EPAD2'],
                                     formats=['>i', dt2dfmt, '>i']))

    outfile = open(outpath, 'wb')
    keys = [key
            for key in ['FLAND', 'VAR1', 'LAI', 'TOPO', 'LUCAT11', 'LUCAT26']
            if key in ncffile.variables.keys()]

    ludts = {'FLAND': _fland_dtype,
             'LUCAT11': _fland_dtype,
             'LUCAT26': _fland_dtype}

    def getempty(key):
        return np.empty(shape=(1,), dtype=ludts.get(key, _other_dtype))

    keyandvar = [(key, getempty(key)) for key in keys]
    for key, var in keyandvar:
        invar = ncffile.variables[key]
        if newstyle:
            if key == 'FLAND'.strip():
                key = 'LUCAT%02d' % nland

            var['SPAD1'] = 8
            var['KEY'] = key.ljust(8)
            var['EPAD1'] = 8

        var['SPAD2'][...] = invar.size * 4
        var['DATA'][...] = invar[:]
        var['EPAD2'][...] = invar.size * 4
        var.tofile(outfile)

    outfile.flush()
    return outfile


registerwriter('camxfiles.landuse', ncf2landuse)
registerwriter('landuse', ncf2landuse)
