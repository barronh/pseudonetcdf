from __future__ import print_function
import numpy as np
from numpy import testing
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.coordutil import gettimes
from datetime import datetime
from collections import OrderedDict
dtype = np.dtype
def _vord(val):
    try:
        vnew = ord(val)
    except TypeError:
        if len(val) == 0:
            vnew = 0.
        else:
            raise
    return vnew

_vordv = np.vectorize(_vord)

def timeit(key):
    from time import time
    if key in timeit.starts:
        started = timeit.starts.pop(key)
        print(key, (time() - started) / 60., 'min')
    else:
        print(key, time(), 'seconds since epoch')
        timeit.starts[key] = time()
timeit.starts = {}

thdtype = dtype([('YYMMDDHHFF', '>10S'), ('LEVEL', '>2S'), ('GRID', '>2S') , ('INDX', '>4S'), ('Z1', '>4S'), ('MB1', '>14S'), ('MB2', '>14S'),
    ('SRCE', '>S4'), ('MGRID', '>S3'), ('VSYS', '>S2'),
    ('POLLAT', '>S7'), ('POLLON', '>S7'), ('REFLAT','>S7'), ('REFLON','>S7'),
    ('GRIDX','>S7'), ('ORIENT','>S7'), ('TANLAT','>S7'),
    ('SYNCHX','>S7'), ('SYNCHY','>S7'), ('SYNCHLAT','>S7'), ('SYNCHLON','>S7'),
    ('RESERVED', '>S7'), ('NX', '>S3'), ('NY', '>S3'), ('NZ', '>S3'), ('VSYS2', '>S2'), ('LENH', '>S4')])
vhdtype = dtype([('YYMMDDHHFF', '>10S'), ('LEVEL', '>2S'), ('grid', '>2S'), ('VKEY', '>4S'), ('EXP', '>4S'), ('PREC', '>14S'), ('VAR1', '>14S')])

stdprops = dict(
    PRSS = ('Pressure at surface', 'hPa'),
    MSLP = ('Pressure at mean sea level', 'hPa'),
    TMPS = ('Temperature at surface', 'K'),
    TPP6 = ('Total precipitation (6-h)', 'm'),
    UMOF = ('U-Momentum flux', 'N/m2'),
    VMOF = ('V-Momentum flux', 'N/m2'),
    SHTF = ('Sfc sensible heat flux', 'W/m2'),
    LTHF = ('Latent heat flux', 'W/m2'),
    DSWF = ('Downward short wave flux', 'W/m2'),
    T02M = ('Temperature at 2 m', 'K'),
    RH2M = ('Relative humidity at 2 m', '%'),
    U10M = ('U-component of wind at 10 m', 'm/s'),
    V10M = ('V-component of wind at 10 m', 'm/s'),
    TCLD = ('Total cloud cover', '%'),
    UWND = ('U wind component (respect to grid)', 'm/s'),
    VWND = ('V wind component (respect to grid)', 'm/s'),
    HGTS = ('Geopotential height', 'gpm*'),
    TEMP = ('Temperature', 'K'),
    WWND = ('Pressure vertical velocity', 'hPa/s'),
    SPHU = ('Specific humidity', 'unknown'),
    )

def readvardef(vheader, out = {}):
    out['keys'] = keys = OrderedDict()
    out['checksums'] = checksums = OrderedDict()
    out['vglvls'] = vglvls = []
    while not np.char.replace(vheader, b' ', b'') == b'':
        vglvl = float(vheader[:6])
        vgvars = int(vheader[6:8])
        vheader = vheader[8:]
        vglvls.append(vglvl)
        keys[vglvl] = []
        for vgvari in range(vgvars):
            vgkey = vheader[:4]
            checksum = int(vheader[4:7])
            keys[vglvl].append(vgkey)
            checksums[vglvl, vgkey] = checksum
            vheader = vheader[8:]
    out['sfckeys'] = keys[vglvls[0]]
    out['laykeys'] = keys[vglvls[1]]
    return out

def writevardef(vglvls, keys, checksums):
    """
    vglvls - iterables of floats [1, .98, ...]
    keys - dictionary of keys by level {vglvli: [vark, ....], ...}
    checksums - integer checsums for each lvl, key -- {(vglvl, vark): checksum, ...}
    """
    vardef = ''
    vgtxts = getvgtxts(vglvls)
    
    for vglvl, vgtxt in zip(vglvls, vgtxts):
        vardef += vgtxt[:6]
        lvlvars = keys[vglvl]
        vardef += '%2d' % len(lvlvars)
        for varkey in lvlvars:
            checksum = checksums[vglvl, varkey]
            vardef += varkey.ljust(4)[:4].decode()
            vardef += '%3d' % checksum
            vardef += ' '
    
    return vardef.ljust(len(vardef) + 108)


def inqarlpackedbit(path):
    f = open(path, 'r')
    fheader = np.fromfile(f, count = 1, dtype = thdtype)[0]
    out = {}
    out['LENH'] = hlen = int(fheader['LENH'])
    out['NX'] = int(fheader['NX'])
    out['NY'] = int(fheader['NY'])
    nz = out['NZ'] = int(fheader['NZ'])
    vheader = np.fromfile(f, count = 1, dtype = '>%dS' % hlen)[0]
    readvardef(vheader, out)
    return out

def getvgtxts(vglvls):
    vgtxts = []
    for vglvl in vglvls:
        if vglvl == 0:
            dp = 0
        else:
            dp = np.floor(np.log10(vglvl)+1)
        tmpl = '%%6.%df' % (5 - dp)
        vgtxt = (tmpl % vglvl)[-6:]
        if len(vgtxt) != 6:
            import pdb; pdb.set_trace()
        vgtxts.append(vgtxt)
    return vgtxts
    
def writearlpackedbit(infile, path):
    """
    path - path to existing arl packed bit file or location for new file
    infile - NetCDF-like file with
        - vertical 2-D (surface) and 3-D (layers) variables 
          with 4 character names
        - a z variable with vertical coordinates
        - all properties from the first index record
    
    """
    requiredkeys = list(thdtype.names)
    for key in requiredkeys:
        getattr(infile, key)
    
    svars = {}
    lvars = {}
    props = {}
    sfckeys = props['sfckeys'] = []
    laykeys = props['laykeys'] = []
    for vark, var in infile.variables.items():
        if len(var.shape) == 3:
            svars[vark] = var
            sfckeys.append(vark.encode('ascii'))
        elif len(var.shape) == 4:
            lvars[vark] = var
            lvar = var
            laykeys.append(vark.encode('ascii'))
        else:
            pass
    
    vglvls = np.append(float(infile.SFCVGLVL), infile.variables['z'][:].array())
    vgtxts = getvgtxts(vglvls)
    # plus one includes surface
    props['NZ'] = lvar.shape[1] + 1
    props['NY'] = lvar.shape[2]
    props['NX'] = lvar.shape[3]
    datamap = maparlpackedbit(path, mode = 'write', shape = (lvar.shape[0],), **props)
    
    theads = datamap['timehead']
    #for key in thdtype.names:
    #    theads[key] = getattr(infile, key)
    
    # vardefs = datamap['vardef']
    # for ti, vardef in enumerate(vardefs):
    #     sfcvardeftxt = vgtxts[0] + '%2d' % len(svars) + ''.join(['%-4s -1 ' % skey.decode() for skey in sfckeys])
    #     layvardeftxt = ''
    #     for vgtxt in vgtxts[1:]:
    #         layvardeftxt += vgtxt + '%2d' % len(lvars) + ''.join(['%-4s -1 ' % lkey.decode() for lkey in laykeys])
    #     
    # 
    # vardeftxt = sfcvardeftxt + layvardeftxt            
    # defsize = 8+8*len(svars)+(8+8*len(lvars))*len(vgtxts[1:])
    # assert(len(vardeftxt) == defsize)
    # vardefs[:] = vardeftxt
    YYMMDDHHFF = getattr(infile, 'YYMMDDHHFF', '0000000000')
    YYMMDDHH = YYMMDDHHFF[:-2]
    FF = YYMMDDHHFF[-2:]
    times = gettimes(infile)
    
    checksums = {}
    for ti, (time, thead) in enumerate(zip(times, theads)):
        for propk in thead.dtype.names:
            thead[propk] = getattr(infile, propk)
        
        timestr = time.strftime('%y%m%d%H').encode('ascii') + FF
        thead['YYMMDDHHFF'] = timestr
        for sfck in sfckeys:
            invar = infile.variables[sfck.decode()]
            var_time = datamap['surface'][sfck.decode()]
            
            
            for varpropk in var_time['head'].dtype.names:
                if not varpropk in ('YYMMDDHHFF', 'LEVEL'):
                    var_time['head'][varpropk][ti] = getattr(invar, varpropk)
            indata = invar[ti]
            
            CVAR, PREC, NEXP, VAR1, KSUM = pack2d(indata, verbose = False) #sfck == b'SNOW')
            var_time['head']['YYMMDDHHFF'][ti] = timestr
            var_time['head']['LEVEL'][ti] = '%2d' % 0
            var_time['head']['PREC'][ti] = '%14.7E' % PREC
            var_time['head']['EXP'][ti] = '%4d' % NEXP
            var_time['head']['VAR1'][ti] = '%14.7E' % VAR1
            checksums[vglvls[0], sfck] = KSUM
            var_time['data'][ti] = CVAR
        for layk in laykeys:
            invar = infile.variables[layk.decode()]
            var_time = datamap['layers'][layk.decode()][ti]
            for li, var_time_lay in enumerate(var_time):
                varhead = var_time_lay['head']            
                for varpropk in varhead.dtype.names:
                    if not varpropk in ('YYMMDDHHFF', 'LEVEL'):
                        varhead[varpropk] = getattr(invar, varpropk)
                varhead['YYMMDDHHFF'] = timestr
                varhead['LEVEL'] = '%2d' % (li + 1)
                indata = invar[ti, li]
                CVAR, PREC, NEXP, VAR1, KSUM = pack2d(indata)
                var_time_lay['data'] = CVAR
                varhead['PREC'] = '%14.7E' % PREC
                varhead['EXP'] = '%4d' % NEXP
                varhead['VAR1'] = '%14.7E' % VAR1
                vglvl = vglvls[li + 1]
                checksums[vglvl, layk] = KSUM
        
        keys = {vglvls[0]: sfckeys}
        for vglvl in vglvls[1:]:
            keys[vglvl] = laykeys
        vardef = writevardef(vglvls, keys, checksums)
        
        datamap['vardef'][ti] = ' '.ljust(datamap['vardef'][ti].itemsize)
        datamap['hdr'][ti] = ' '.ljust(datamap['hdr'][ti].itemsize)
        datamap['vardef'][ti] = vardef.encode('ascii')
    
    
        
    
def maparlpackedbit(path, mode = 'r', shape = None, **props):
    """
    pg 4-9 of http://niwc.noaa.inel.gov/EOCFSV/hysplit/hysplituserguide.pdf
    
    path - path to existing arl packed bit file or location for new file
    props - if props is empty, the arl path must exist. Otherwise:
        required:
            - NZ integer number of layers including surface
            - NY integer number of rows
            - NX integer number of cols
            - sfckeys variable names in the surface layer
            - laykeys variable names in the layers above the surface
        optional:
            - vglvls list nz of vertical coordinates level
            - checksums dictionary {(vglvl, varkey): integer checksum, ...}
    
    FILE FORMAT:
        recl = 50 + NX*NY
        for time in times:
            1 rec of length recl with meta-data (thdtype + vardefdtype + hdrdtype)
            for sfckey in sfckeys:
                1 rec of length recl (vhdtype + nx*ny bytes)
            for layer in layers:
                for laykey in laykeys:
                    1 rec of length recl (vhdtype + nx*ny bytes)
    """
    #import pdb; pdb.set_trace()
    if props == {}:
        props = inqarlpackedbit(path)
    else:
        srflen = 6 + 2 + (4 + 3 + 1) * len(props['sfckeys'])
        laylen = (6 + 2 + (4 + 3 + 1) * len(props['laykeys'])) * (props['NZ'] - 1)
        props['LENH'] = 108 + srflen + laylen
    
    nx = props['NX']
    ny = props['NY']
    # minus 1 excludes surface
    nz = props['NZ'] - 1
    hlen = props['LENH']
    sfckeys = props['sfckeys']
    laykeys = props['laykeys']
    ncell = nx * ny
    vardefdtype = dtype('%d>S' % hlen)
    hdrdtype = dtype('%d>c' % (50 + ncell - hlen - thdtype.itemsize))
    lay1dtype = dtype([('head', vhdtype), ('data', dtype('(%d,%d)>1S' % (ny, nx)))])
    sfcdtype = dtype(dict(names = [k.decode() for k in sfckeys], formats = [lay1dtype] * len(sfckeys)))
    laydtype = dtype(dict(names = [k.decode() for k in laykeys], formats = [lay1dtype] * len(laykeys)))
    layersdtype = dtype([('layers', laydtype, (nz, ))])
    timedtype = dtype([('timehead', thdtype), ('vardef', vardefdtype), ('hdr', hdrdtype), ('surface', sfcdtype), ('layers', laydtype, (nz, ))])
    datamap = np.memmap(path, timedtype, shape = shape, mode = mode)
    return datamap

def unpack(bytes, VAR1, EXP):
    """
    bytes - nx by ny bytes
    VAR1 - as read directly from LABEL as BYTES with dimension time 
    EXP - EXP as read directly from LABEL as BYTES with dimension time 
    """
    vold = VAR1.astype('f')
    scale = np.float32(2.0)**np.float32(7-EXP.astype('i'))
    invscale = np.float32(1.) / scale
    data = (bytes.view('uint8') - np.float32(127.)) * invscale[..., None, None]
    data[..., 0, 0] += vold
    data[..., 0] = np.cumsum(data[..., 0], axis = data.ndim - 2)
    vdata = np.cumsum(data, axis = data.ndim - 1)
    return vdata

def CHAR(c):
    return chr(c).encode('ascii')
    

def pack2d(RVARA, verbose = False):
    """
    CHARACTER, INTENT(OUT) :: cvar(nxy)   ! packed char*1 output array
    REAL,      INTENT(OUT) :: prec        ! precision of packed data array
    INTEGER,   INTENT(OUT) :: nexp        ! packing scaling exponent
    REAL,      INTENT(OUT) :: var1        ! value of real array at position (1,1)
    INTEGER,   INTENT(OUT) :: ksum        ! rotating checksum of packed data 
    """
    MAX = np.maximum
    FLOAT = np.float32
    INT = np.int32
    ABS = np.abs
    LOG = np.log
    NY, NX = RVARA.shape
    NXY = NX * NY
    CVAR = np.zeros(RVARA.shape, dtype = 'uint8')
    RVAR = RVARA.astype('f')
    VAR1=RVAR[0,0]
    
    # find the maximum difference between adjacent elements
    # START ORIGINAL SERIAL CODE
    # ROLD= VAR1
    # RMAX= 0.0
    #for J in range(NY):
    #    for I in range(NX):
    #        # compute max difference between elements along row
    #        RMAX=MAX( ABS(RVAR[J, I]-ROLD), RMAX)
    #        ROLD=RVAR[J,I]
    #    
    #    # row element 1 difference always from previous row
    #    ROLD=RVAR[J, 0]
    #
    # END ORIGINAL SERIAL CODE
    # START NUMPY VECTOR CODE
    colmax = np.abs(np.diff(RVAR, axis = 1)).max()
    rowmax = np.abs(np.diff(np.append(VAR1, RVAR[:, 0]), axis = 0)).max()
    RMAX = np.maximum(colmax, rowmax)
    # END NUMPY VECTOR CODE
    
    SEXP=0.0
    # compute the required scaling exponent
    if RMAX != 0.0:
        SEXP=LOG(RMAX)/LOG(np.float32(2.))
    
    NEXP=INT(SEXP)
    # positive or whole number scaling round up for lower precision
    if SEXP >= 0.0 or (SEXP % 1.0) == 0.0:
        NEXP=NEXP+1
    # precision range is -127 to 127 or 254
    PREC=np.float32((2.0**NEXP)/254.0)
    SCEXP=np.float32(2.0**(7-NEXP))
    
    # START ORIGINAL SERIAL CODE
    # # initialize checksum
    # KSUM=0
    # 
    # K=0
    # # set column1 value
    # RCOL=VAR1
    # RCOL = VAR1
    # 
    # # pack the array from low to high
    # for J in range(NY):
    #     ROLD = RCOL
    #     for I in range(NX):
    #         K=K+1
    #         # packed integer at element
    #         ICVAL=INT((RVAR[J, I]-ROLD)*SCEXP+127.5)
    #         
    #         # convert to character
    #         #CVAR[J, I]=CHAR(ICVAL)
    #         CVAR[J, I]=ICVAL
    #         # previous element as it would appear unpacked
    #         ROLD=FLOAT(ICVAL-127)/SCEXP+ROLD
    #         # save the first column element for next row
    #         if I == 0:
    #             RCOL=ROLD.copy()
    #         # maintain rotating checksum
    #         KSUM=KSUM+ICVAL
    #         # if sum carries over the eighth bit add one
    #         if KSUM >= 256:
    #             KSUM=KSUM-255
    # 
    # CVART = CVAR.copy()
    # KSUMT = KSUM
    # END ORIGINAL SERIAL CODE

    # START NUMPY VECTOR CODE
    ROLDS = np.zeros_like(RVAR[:, 0])
    ROLD = VAR1
    for J in range(NY):
        ICVAL = INT((RVAR[J, 0] - ROLD) * SCEXP + 127.5)
        CVAR[J, 0] = ICVAL
        ROLD=FLOAT(ICVAL-127)/SCEXP+ROLD
        ROLDS[J] = ROLD
    
    ROLD = ROLDS
    for I in range(1, NX):
        ICVAL = INT((RVAR[:, I] - ROLD)*SCEXP+127.5)
        CVAR[:, I] = ICVAL
        ROLD = FLOAT(ICVAL - 127)/SCEXP+ROLD
    KSUM = INT(CVAR.sum()) % 255
    # END NUMPY VECTOR CODE
    
    #assert((CVART == CVAR).all())
    #assert(KSUM == KSUMT)
    return CVAR.view('>S1'), PREC, NEXP, VAR1, KSUM

class arlpackedbit(PseudoNetCDFFile):
    """
    Format as follows:
    
    for t in times:
        thdtype = dtype([('YYMMDDHHFF', '>10S'), ('LEVEL', '>2S'), ('GRID', '>2S') , ('INDX', '>4S'), ('Z1', '>4S'), ('MB', '2>14S'),  ('SRCE', '>S4'), ('MGRID', '>S3'), ('VSYS', '>S2'), ('POLLAT', '>S7'), ('POLLON', '>S7'), ('REFLAT','>S7'), ('REFLON','>S7'), ('GRIDX','>S7'), ('ORIENT','>S7'), ('TANLAT','>S7'), ('SYNCHX','>S7'), ('SYNCHY','>S7'), ('SYNCHLAT','>S7'), ('SYNCHLON','>S7'), ('RESERVED', '>S7'), ('NX', '>S3'), ('NY', '>S3'), ('NZ', '>S3'), ('VSYS2', '>S2'), ('LENH', '>S4')]) + ny * nx - 158 - 50
        ny+nx-
        for surfacekey in surfacekeys:
            YYMMDDHHFFLLGGVVVVEEEEPPPPPPPPPPPPPPVVVVVVVVVVVVVV+
            P(ny,nx)
        for layerkey in layerkeys:
            for layer in layers:
                YYMMDDHHFFLLGGVVVVEEEEPPPPPPPPPPPPPPVVVVVVVVVVVVVV+
                P(ny,nx)
                
    
    YY = %y in strftime for GMT date
    MM = %m in strftime
    DD = %d in strftime
    HH = %H of model hour
    FF = %H of forecast hour
    LL = 2-digit layer
    GG = 2-digit grid
    IIII = 'INDX'
    ZZZZ = 4-digit integer 0
    XXXXXXXXXXXXX = %14.7e 0
    YYYYYYYYYYYYY = %14.7e 0
    EEEE = 4-digit integer exponent
    PPPPPPPPPPPPP = %14.7e precision
    VVVVVVVVVVVVV = %14.7e value R(1,1)
    P(ny,nx) = ny by nx 1-byte values dervied from real values (R) according to the formula below.
    
    P(i,j) = (Ri,j  - Ri-1,j)* (2**(7-(ln dRmax / ln 2)))
    """
    def __init__(self, path):
        self._path = path
        self._f = f = open(path, 'r') 
        f.seek(0, 2) 
        fsize = f.tell()
        f.seek(0, 0) 

        vord = np.vectorize(ord)
#        timesfc = []
#        times = []
#        tflag = []
        
        f.seek(0, 0) 
        datamap = self._datamap = maparlpackedbit(path)
        t0hdr = datamap['timehead'][0]
        for key in thdtype.names:
            setattr(self, key, t0hdr[key])

        out = readvardef(datamap['vardef'][0])
        self.XSIZE = self.YSIZE = float(t0hdr['GRIDX']) * 1000.
        self.createDimension('time', datamap.shape[0])
        self.createDimension('x', int(self.NX))
        self.createDimension('y', int(self.NY))
        # minus 1 excludes surface
        self.createDimension('z', int(self.NZ) - 1)

        tflag = np.char.replace(datamap['timehead']['YYMMDDHHFF'], b' ', b'0')
        sfckeys = self._datamap['surface'].dtype.names
        laykeys = self._datamap['layers'][0].dtype.names
        
        self.variables = PseudoNetCDFVariables(func = self._getvar, keys = list(sfckeys + laykeys))
        times = [datetime.strptime(t.astype('S8').decode(), '%y%m%d%H') for t in tflag]
        rtime = times[0]
        hours_since = [(t - rtime).total_seconds() // 3600 for t in times]
        timev = self.createVariable('time', 'i', ('time',))
        timev[:] = hours_since
        timev.units = rtime.strftime('hours since %F %H:%M:%S')
        z = self.createVariable('z', 'f', ('z',))
        z[:] = out['vglvls'][1:]
        self.SFCVGLVL = out['vglvls'][0]
        z.units = {1: 'pressure sigma', 2: 'pressure absolute', 3: 'terrain sigma', 4: 'hybrid sigma'}.get(int(self.VSYS2), 'unknown')
        x = self.createVariable('x', 'f', ('x',))
        x[:] = np.arange(x[:].size) * float(t0hdr['GRIDX']) * 1000.
        y = self.createVariable('y', 'f', ('y',))
        y[:] = np.arange(y[:].size) * float(t0hdr['GRIDX']) * 1000.

    
    def _getvar(self, k):
        datamap = self._datamap
        stdname, stdunit = stdprops.get(k, (k, 'unknown'))
        if k in datamap['surface'].dtype.names:
            vhead = datamap['surface'][k]['head']
            v11 = vhead['VAR1']
            EXP = vhead['EXP']
            props = dict([(k, vhead[k][0]) for k in vhead.dtype.names if not k in ('YYMMDDHHFF', 'LEVEL')])
            bytes = datamap['surface'][k]['data']
            vdata = unpack(bytes, v11, EXP)
            #if k == 'MXLR':
            #    CVAR, PREC, NEXP, VAR1, KSUM = pack2d(vdata[21], verbose = True)
                
                    
            return PseudoNetCDFVariable(self, k, 'f', ('time', 'y', 'x'), values = vdata, units = stdunit, standard_name = stdname, **props)
        elif k in datamap['layers'].dtype.names:
            bytes = datamap['layers'][k]['data']
            vhead = datamap['layers'][k]['head']
            v11 = vhead['VAR1']
            EXP = vhead['EXP']
            props = dict([(k, vhead[k][0, 0]) for k in vhead.dtype.names if not k in ('YYMMDDHHFF', 'LEVEL')])
            vdata = unpack(bytes, v11, EXP)
            return PseudoNetCDFVariable(self, k, 'f', ('time', 'z', 'y', 'x'), values = vdata, units = stdunit, standard_name = stdname, **props)

from PseudoNetCDF._getwriter import registerwriter
registerwriter('noaafiles.arlpackedbit', writearlpackedbit)
registerwriter('arlpackedbit', writearlpackedbit)

if __name__ == '__main__':
    import sys
    out = arlpackedbit(sys.argv[1])
