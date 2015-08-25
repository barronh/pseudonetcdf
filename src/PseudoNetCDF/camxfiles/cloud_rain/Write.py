import numpy as np
def ncf2cloud_rain(ncffile, outpath, tflag = 'TFLAG'):
    outfile=open(outpath,'wb')
    varkeys = [vk for vk in ['CLOUD', 'PRECIP', 'RAIN', 'SNOW', 'GRAUPEL', 'COD'] if vk in ncffile.variables.keys()]
    buf = np.array(len(ncffile.FILEDESC) + 12, dtype = '>i').tostring()
    outfile.write(buf)
    outfile.write(np.array(ncffile.FILEDESC, dtype = 'c'))
    nxcl = len(ncffile.dimensions['COL'])
    nycl = len(ncffile.dimensions['ROW'])
    nzcl = len(ncffile.dimensions['LAY'])
    outfile.write(np.array([nxcl, nycl, nzcl], dtype = '>i'))
    outfile.write(buf)
        
    for di, (d,t) in enumerate(ncffile.variables[tflag][:,0,:]):
        t=np.array(t.astype('>f')/100,ndmin=1).astype('>f')
        d=np.array(d,ndmin=1).astype('>i')
        d=(d%(d/100000*100000)).astype('>i')
        buf = np.array(8, dtype = '>i').tostring()
        outfile.write(buf+t.tostring()+d.tostring()+buf)
        for zi in range(nzcl):
            for varkey in varkeys:
                vals = ncffile.variables[varkey][di, zi].astype('>f')
                buf=np.array((vals.size)*4,ndmin=1).astype('>i').tostring()
                outfile.write(buf)
                vals.tofile(outfile)
                outfile.write(buf)
    
    outfile.flush()
    return outfile
