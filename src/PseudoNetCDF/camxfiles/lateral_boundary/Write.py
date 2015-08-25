import numpy as np
_emiss_hdr_fmt=np.dtype(dict(names=['SPAD','name','note','itzon','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','(10,4)>S1','(60,4)>S1','>i','>i','>i','>f','>i','>f','>i']))

_grid_hdr_fmt=np.dtype(dict(names=['SPAD','plon','plat','iutm','xorg','yorg','delx','dely','nx','ny','nz','iproj','istag','tlat1','tlat2','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))

_cell_hdr_fmt=np.dtype(dict(names=['SPAD','ione1','ione2','nx','ny','EPAD'],formats=['>i','>i','>i','>i','>i','>i']))

_time_hdr_fmt=np.dtype(dict(names=['SPAD','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))

_spc_fmt=np.dtype("(10,4)>S1")

def ncf2lateral_boundary(ncffile, outpath):
    emiss_hdr = np.zeros(shape = (1,), dtype = _emiss_hdr_fmt)
    emiss_hdr[0]['name'][:, :] = ' '
    emiss_hdr[0]['name'][:, 0] = np.array(ncffile.NAME, dtype = '>c')
    emiss_hdr[0]['note'][:, :] = ' '
    emiss_hdr[0]['note'][:, 0] = np.array(ncffile.NOTE, dtype = '>c')
    gdtype = getattr(ncffile, 'GDTYPE', -999)
    emiss_hdr['itzon'][0] = ncffile.ITZON
    nspec = len(ncffile.dimensions['VAR']) / 4
    emiss_hdr['nspec'] = nspec
    emiss_hdr['ibdate'] = ncffile.SDATE%(ncffile.SDATE/100000*100000)
    emiss_hdr['btime'] = ncffile.STIME / 100.
    EDATE, ETIME = ncffile.variables['TFLAG'][-1,0,:]
    emiss_hdr['iedate'] = EDATE%(EDATE/100000*100000)
    emiss_hdr['etime'] = ETIME / 10000. + 1.
    emiss_hdr['SPAD'] = _emiss_hdr_fmt.itemsize - 8
    emiss_hdr['EPAD'] = _emiss_hdr_fmt.itemsize - 8
    
    NCOLS = len(ncffile.dimensions['COL'])
    NROWS = len(ncffile.dimensions['ROW'])
    NLAYS = len(ncffile.dimensions['LAY'])
    grid_hdr = np.zeros(shape = (1,), dtype = _grid_hdr_fmt)
    grid_hdr['SPAD'] = grid_hdr.itemsize - 8
    grid_hdr['plon'] = ncffile.PLON
    grid_hdr['plat'] = ncffile.PLAT
    grid_hdr['iutm'][0] = ncffile.IUTM
    grid_hdr['xorg'] = ncffile.XORIG
    grid_hdr['yorg'] = ncffile.YORIG
    grid_hdr['delx'] = ncffile.XCELL
    grid_hdr['dely'] = ncffile.YCELL
    grid_hdr['nx'] = NCOLS
    grid_hdr['ny'] = NROWS
    grid_hdr['nz'] = NLAYS
    grid_hdr['iproj'] = ncffile.CPROJ
    grid_hdr['tlat1'] = ncffile.TLAT1
    grid_hdr['tlat2'] = ncffile.TLAT2
    grid_hdr['istag'] = ncffile.ISTAG
    grid_hdr['rdum5'] = 0.
    grid_hdr['EPAD'] = grid_hdr.itemsize - 8

    cell_hdr = np.zeros(shape = (1,), dtype = _cell_hdr_fmt)
    cell_hdr['SPAD'] = cell_hdr['EPAD'] = cell_hdr.itemsize - 8
    cell_hdr['ione1'] = 1
    cell_hdr['ione2'] = 1
    cell_hdr['nx'] = grid_hdr['nx']
    cell_hdr['ny'] = grid_hdr['ny']

    time_hdr = np.zeros(shape = (len(ncffile.dimensions['TSTEP']),), dtype = _time_hdr_fmt)
    time_hdr['SPAD'] = 16
    time_hdr['EPAD'] = 16
    date, time = ncffile.variables['TFLAG'][:, 0].T
    time = time.astype('>f') / 10000.
    date = date%(date/100000*100000)
    time_hdr['ibdate'] = date
    time_hdr['btime'] = time
    time_hdr['iedate'] = date
    time_hdr['etime'] = time + 1.
    time_hdr['iedate'] += time_hdr['etime'] // 24
    time_hdr['etime'] -= (time_hdr['etime'] // 24) * 24
    emiss_hdr['ibdate'] = time_hdr['ibdate'][0]
    emiss_hdr['btime'] = time_hdr['btime'][0]
    emiss_hdr['iedate'] = time_hdr['iedate'][-1]
    emiss_hdr['etime'] = time_hdr['etime'][-1]
    emiss_hdr['SPAD'] = _emiss_hdr_fmt.itemsize - 8
    emiss_hdr['EPAD'] = _emiss_hdr_fmt.itemsize - 8
    
    
    spc_hdr = np.zeros(shape = (1,), dtype = dict(names = ['SPAD1', 'DATA', 'EPAD1'], formats = ['>i', np.dtype("(%d,10,4)>S1" % nspec), '>i']))
    spc_hdr['SPAD1'] = nspec * 40
    spc_hdr['EPAD1'] = nspec * 40
    spc_names = np.array(getattr(ncffile, 'VAR-LIST'), dtype = '>c').reshape(-1, 16)[:-1:4, 5:-1].copy()
    spc_hdr[0]['DATA'][:] = ' '
    spc_hdr[0]['DATA'][:, :, 0] = spc_names
    spc_names = spc_names.view('>S10')
    nz = len(ncffile.dimensions['LAY'])
    outfile = open(outpath, 'wb')
    emiss_hdr.tofile(outfile)
    grid_hdr.tofile(outfile)
    cell_hdr.tofile(outfile)
    spc_hdr.tofile(outfile)
    for ename, ei in [('WEST', 1), ('EAST', 2), ('SOUTH', 3), ('NORTH', 4)]:
        if hasattr(ncffile, '_boundary_def'):
            ncffile._boundary_def[ename].tofile(outfile)
        else:
            nbcell = dict(WEST = NROWS, EAST = NROWS,
                          SOUTH = NCOLS, NORTH = NCOLS)[ename]
            buf = (nbcell * 4 + 3) * 4
            np.array([buf, 1, ei, nbcell, 0, 0, 0, 0] + [2, 0, 0, 0] * (nbcell - 2) + [0, 0, 0, 0, buf]).astype('>i').tofile(outfile)
    
    for di, (d, t) in enumerate(ncffile.variables['TFLAG'][:, 0]):
        time_hdr[di].tofile(outfile)
        for spc_key, spc_name in zip(spc_names, spc_hdr[0]['DATA']):
            for ename, ei in [('WEST', 1), ('EAST', 2), ('SOUTH', 3), ('NORTH', 4)]:
                var = ncffile.variables[ename + '_' + str(np.char.strip(spc_key[0]))]
                data = var[di].astype('>f')
                buf = np.array(4+40+4+data.size*4).astype('>i')
                buf.tofile(outfile)
                np.array(1).astype('>i').tofile(outfile)
                spc_name.tofile(outfile)
                np.array(ei).astype('>i').tofile(outfile)
                data.tofile(outfile)
                buf.tofile(outfile)
    
    outfile.flush()
    return outfile
