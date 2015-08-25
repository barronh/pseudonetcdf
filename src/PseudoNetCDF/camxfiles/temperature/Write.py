import numpy as np
"""
hour,idate,((temps(i,j),i=1,nx),j=1,ny)
Loop from 1 to nlay layers:
    hour,idate,((temp(i,j,k),i=1,nx),j=1,ny)
"""

def ncf2temperature(ncffile, outpath):
    outfile = open(outpath, 'wb')
    sfc = ncffile.variables['SURFTEMP']
    air = ncffile.variables['AIRTEMP']
    nz, nr, nc = air.shape[-3:]
    nelem = nr * nc * 4 + 8
    for di,(d,t) in enumerate(ncffile.variables['TFLAG'][:, 0]):
        t=np.array(t/100,ndmin=1, dtype = '>f')
        d=np.array(d,ndmin=1).astype('>i')
        d=(d%(d/100000*100000)).astype('>i')
        buf = np.array(nelem, dtype = '>i').tostring()
        outfile.write(buf)
        t.tofile(outfile)
        d.tofile(outfile)
        sfc[di].astype('>f').tofile(outfile)
        outfile.write(buf)
        for zi in range(nz):
            outfile.write(buf)
            t.tofile(outfile)
            d.tofile(outfile)
            air[di, zi].astype('>f').tofile(outfile)
            outfile.write(buf)
    
    outfile.flush()
    return outfile
