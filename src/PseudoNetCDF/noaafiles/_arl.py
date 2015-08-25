from numpy import *
from numpy import testing
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables

def _vord(val):
    try:
        vnew = ord(val)
    except TypeError:
        if len(val) == 0:
            vnew = 0.
        else:
            raise
    return vnew

_vordv = vectorize(_vord)

def timeit(key):
    from time import time
    if key in timeit.starts:
        started = timeit.starts.pop(key)
        print key, (time() - started) / 60., 'min'
    else:
        print key, time(), 'seconds since epoch'
        timeit.starts[key] = time()
timeit.starts = {}

thdtype = dtype([('YYMMDDHHMM', '>10S'), ('LEVEL', '>2S'), ('grid', '>2S') , ('INDX', '>4S'), ('z1', '>4S'), ('mb', '2>14S'), ('SRCE', '>S4'), ('MGRID', '>S3'), ('VSYS', '>S2'), ('POL', '(2,)>S7'), ('REF','(2,)>S7'), ('GRIDX','>S7'), ('ORIENT','>S7'), ('LAT0','>S7'), ('EQUATE','(4,)>S7'), ('dum', '>S7'), ('nx', '>S3'), ('ny', '>S3'), ('nz', '>S3'), ('VSYS2', '>S2'), ('LENH', '>S4')])
vhdtype = dtype([('YYMMDDHH', '>10S'), ('LEVEL', '>2S'), ('grid', '>2S'), ('VKEY', '>4S'), ('EXP', '>4S'), ('PREC', '>14S'), ('VAR1', '>14S')])

class arlpackedbit(PseudoNetCDFFile):
    def _readindx(self):
        f = self._f
        self._fheader = fheader = fromfile(f, count = 1, dtype = thdtype)[0]
        hlen = int(fheader['LENH'])
        vheader = fromfile(f, count = 1, dtype = '>%dS' % hlen)[0]
        self._keys = keys = {}
        self._vglvls = vglvls = []
        nx = self._nx = int(fheader['nx'])
        ny = int(fheader['ny'])
        ncell = nx*ny
        self.NROWS = ny
        self.NCOLS = nx
        self.NLAYS = int(fheader['nz']) - 1
        while len(vheader) > 0:
            vglvl = float(vheader[:6])
            vgvars = int(vheader[6:8])
            vheader = vheader[8:]
            vglvls.append(vglvl)
            keys[vglvl] = []
            for vgvari in range(vgvars):
                vgkey = vheader[:4]
                checksums = int(vheader[4:7])
                keys[vglvl].append((vgkey, checksums))
                vheader = vheader[8:]

        data = fromfile(f, count = (50 + ncell - hlen - fheader.dtype.itemsize), dtype = 'c')
        self._varrec = dtype('%d>S' % hlen)
        self._hdrrec = dtype('%d>c' % (50 + ncell - hlen - fheader.dtype.itemsize))
        

    def __init__(self, path):
        self._f = f = open(path, 'r') 
        f.seek(0, 2) 
        fsize = f.tell()
        f.seek(0, 0) 

        vord = vectorize(ord)
#        timesfc = []
#        times = []
#        tflag = []
        self._readindx()
        fheader = self._fheader
        nz, ny, nx = self.NLAYS, self.NROWS, self.NCOLS
        f.seek(0, 0) 
        self._sfckeys = sfckeys = [k for k, c in self._keys[self._vglvls[0]]]
        self._laykeys = laykeys = [k for k, c in self._keys[self._vglvls[1]]]
        lay1dtype = dtype([('head', vhdtype), ('data', dtype('(%d,%d)>1S' % (ny, nx)))])
        sfcdtype = dtype(dict(names = sfckeys, formats = [lay1dtype] * len(sfckeys)))
        laydtype = dtype(dict(names = laykeys, formats = [lay1dtype] * len(laykeys)))
        layersdtype = dtype([('layers', laydtype, (nz, ))])
        timedtype = dtype([('timehead', thdtype), ('varrec', self._varrec), ('hdrrec', self._hdrrec), ('surface', sfcdtype), ('upper', layersdtype)])
        self._datamap = datamap = memmap(path, timedtype, mode = 'r')
#         timeit('memmap all')
#         timeit('fromfile t')
#             
#         # Outer loop is time
#         while f.tell() != fsize:
#             # Read the time index
#             self._readindx()
#             fheader = self._fheader
#             ny, nx = self.NROWS, self.NCOLS
#             vdata = zeros((ny, nx), dtype = 'f')
#             YYMMDDHHMM = fheader[0].replace(' ', '0')
#             tflag.append(YYMMDDHHMM)
#             timesfc.append({})
#             times.append({})
#             for li, levelkey in enumerate(self._vglvls):
#                 for key, chksum in self._keys[levelkey]:
#                     vdata[:] = 0.
#                     if f.tell() == fsize: break
#                     vheader = fromfile(f, count = 1, dtype = vhdtype)[0]
#                     scale = 2.0**(7-int(vheader['EXP']))
#                     invscale = 1./scale
#                         
# 
#                     #print vheader['YYMMDDHH'], vheader['LEVEL'], key
#                     assert(li == int(vheader['LEVEL']))
#                     data = fromfile(f, count = 1, dtype = '(%d,%d)>1S' % (ny, nx))[0]
#                     vdatat = ((_vordv(data) - 127.) * invscale).astype('f')
#                     vdatat[0, 0] = float(vheader['VAR1'])
#                     vdatat[:, 0] = cumsum(vdatat[:, 0])
#                     vdatat = cumsum(vdatat, axis = 1)
#                     #for ri, row in enumerate(data):
#                     #    if ri == 0:
#                     #        vold = float(vheader['VAR1'])
#                     #    else:
#                     #        vold = vdata[ri-1, 0]
#                     #    rowvals = (_vordv(row) - 127.) * invscale
#                     #    rowvals[0] += vold
#                     #    rowvalsc = cumsum(rowvals).astype('f')
#                     #    vdata[ri, :] = rowvalsc
#                     #testing.assert_almost_equal(vdatat, vdata, decimal=4, err_msg='Ahh', verbose=True)
#                     vdata = vdatat
#                     if li == 0:
#                         timesfc[-1][key] = vdata.copy()
#                     else:
#                         if key not in times[-1]:
#                             times[-1][key] = vdata[None,:, :]
#                         else:
#                             times[-1][key] =  concatenate([times[-1][key], vdata[None,:, :]], axis = 0)
#             break
#         timeit('fromfile t')
#         for k in timesfc[0].keys():
#             outvar = self.createVariable(k, 'f', ('time', 'y', 'x'), values = array([d[k] for d in timesfc]))
#             outvar.long_name = k
# 
#         for k in times[0].keys():
#             outvar = self.createVariable(k, 'f', ('time', 'z', 'y', 'x'), values = array([d[k] for d in times]))
#             outvar.long_name = k
        self.YCENT, self.XCENT = map(float, datamap['timehead']['POL'][0])
        self.YCENT_ROW, self.XCENT_COL, self.YCENT_LAT, self.XCENT_LON = map(float, fheader[14])
        self.XSIZE = self.YSIZE = float(datamap['timehead']['GRIDX'][0]) * 1000.
        self.createDimension('x', nx)
        self.createDimension('y', ny)
        self.createDimension('z', self.NLAYS)
        self.createDimension('time', datamap.shape[0])
        tflag = char.replace(datamap['timehead']['YYMMDDHHMM'], ' ', '0')
        self.variables = PseudoNetCDFVariables(func = self._getvar, keys = sfckeys + laykeys)
        self.createVariable('TFLAG', 'S10', ('time',), values = array(tflag, dtype = 'S10'))
    def _getvar(self, k):
        datamap = self._datamap
        if k in self._sfckeys:
            vold = datamap['surface'][k]['head']['VAR1'].astype('f')
            scale = 2.0**(7-datamap['surface'][k]['head']['EXP'].astype('i'))
            invscale = 1. / scale
            data = (datamap['surface'][k]['data'].view('uint8') - 127.) * invscale[:, None, None]
            data[:, 0, 0] += vold
            data[:, :, 0] = cumsum(data[:, :, 0], axis = 1)
            vdata = cumsum(data, axis = 2)
            return PseudoNetCDFVariable(self, k, 'f', ('time', 'y', 'x'), values = vdata, units = 'unknown')
        elif k in self._laykeys:
            vold = datamap['upper']['layers'][k]['head']['VAR1'].astype('f')
            scale = 2.0**(7-datamap['upper']['layers'][k]['head']['EXP'].astype('i'))
            invscale = 1. / scale
            data = (datamap['upper']['layers'][k]['data'].view('uint8') - 127.) * invscale[:, :, None, None]
            data[:, :, 0, 0] += vold
            data[:, :, :, 0] = cumsum(data[:, :, :, 0], axis = 2)
            vdata = cumsum(data, axis = 3)
            return PseudoNetCDFVariable(self, k, 'f', ('time', 'z', 'y', 'x'), values = vdata, units = 'unknown')
        
if __name__ == '__main__':
    import sys
    out = arlpackedbit(sys.argv[1])
    print out.XCENT, out.XCENT_COL
    print out.YCENT, out.YCENT_ROW