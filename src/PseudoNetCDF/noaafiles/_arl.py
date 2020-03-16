from __future__ import print_function
import numpy as np
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable
from PseudoNetCDF import PseudoNetCDFVariables
from PseudoNetCDF.coordutil import gettimes
from PseudoNetCDF._getwriter import registerwriter
from datetime import datetime
from collections import OrderedDict
from warnings import warn
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

thdtype = dtype([('YYMMDDHHFF', '>10S'), ('LEVEL', '>2S'), ('GRID', '>2S'),
                 ('INDX', '>4S'), ('Z1', '>4S'), ('MB1', '>14S'),
                 ('MB2', '>14S'), ('SRCE', '>S4'), ('MGRID', '>S3'),
                 ('VSYS', '>S2'), ('POLLAT', '>S7'), ('POLLON', '>S7'),
                 ('REFLAT', '>S7'), ('REFLON', '>S7'),
                 ('GRIDX', '>S7'), ('ORIENT', '>S7'), ('TANLAT', '>S7'),
                 ('SYNCHX', '>S7'), ('SYNCHY', '>S7'), ('SYNCHLAT', '>S7'),
                 ('SYNCHLON', '>S7'), ('RESERVED', '>S7'), ('NX', '>S3'),
                 ('NY', '>S3'), ('NZ', '>S3'), ('VSYS2', '>S2'),
                 ('LENH', '>S4')])
vhdtype = dtype([('YYMMDDHHFF', '>10S'), ('LEVEL', '>2S'), ('grid', '>2S'),
                 ('VKEY', '>4S'), ('EXP', '>4S'), ('PREC', '>14S'),
                 ('VAR1', '>14S')])

stdprops = dict(
    x=('x', 'm'),
    y=('y', 'm'),
    longitude=('longitude', 'degrees'),
    latitude=('latitude', 'degrees'),
    longitude_bounds=('longitude_bounds', 'degrees'),
    latitude_bounds=('latitude_bounds', 'degrees'),
    ABSV=('ABSOLUTE VORTICITY', '100000.00 ? /s'),
    CAPE=('CONVECTIVE AVAILABLE POTENTIAL ENERGY', 'J/kg'),
    CFZR=('CATEGORIAL FREEZING RAIN (YES=1/NO=0)', '0--1'),
    CICE=('CATEGORIAL ICE PELLETS (YES=1/NO=0)', '0--1'),
    CINH=('CONVECTIVE INHIBITION', 'J/kg'),
    CLDB=('CLOUD BOTTOM HEIGHT', 'm'),
    CLDT=('CLOUD TOP HEIGHT', 'm'),
    CONC=('CONCENTRATION', '/m3'),
    CPP3=('3-HOUR ACC. CONVECTIVE PRECIPITATION', 'mm'),
    CPP6=('6-HOUR ACC. CONVECTIVE PRECIPITATION', 'mm'),
    CPPT=('CONVECTIVE PRECIPITATION', 'mm'),
    CRAI=('CATEGORIAL RAIN (YES=1/NO=0)', '0--1'),
    CSNO=('CATEGORIAL SNOW (YES=1/NO=0)', '0--1'),
    CWMR=('CLOUD WATER MIXING RATIO', 'g/kg'),
    DEPO=('SURFACE DEPOSITION', '/m2'),
    DEPS=('SURFACE DEPOSITION', '/m2'),
    DEPV=('DEPOSITION VELOCITY', 'm/s'),
    DLWF=('DOWNWARD LONG WAVE RADIATION FLUX', 'W/m2'),
    DP2M=('2 M DEW POINT', 'degC'),
    DSWF=('DOWNWARD SHORT WAVE RADIATION FLUX', 'W/m2'),
    EXCO=('EXCHANGE COEFFICIENT', '1000.00 ? GM2S'),
    FLAG=('WIND FLAGS', 'KNTS'),
    HCLD=('HIGH CLOUD COVER', 'PCT'),
    HFLX=('LATENT HEAT FLUX', 'W/m2'),
    HGT1=('1000 MB HEIGHT', 'm'),
    HGT5=('500 MB HEIGHT', 'm'),
    HGTS=('Geopotential height', 'gpm*'),
    ICWT=('ICE-COVERED WATER', '0--1'),
    LCLD=('LOW CLOUD COVER', 'PCT'),
    LHTF=('LATENT HEAT NET FLUX', 'W/m2'),
    LIB4=('BEST 4-LAYER LIFTED INDEX', 'K'),
    LISD=('STANDARD LIFTED INDEX', 'degC'),
    LTHF=('Latent heat flux', 'W/m2'),
    LTNG=('LIGHTNING STRIKES', '#'),
    MCLD=('MEDIUM CLOUD COVER', 'PCT'),
    MOMF=('TOTAL MOMENTUM FLUX', 'N/m2'),
    MSLP=('MEAN SEA-LEVEL PRESSURE', 'hPa'),
    MXHT=('MIXED LAYER DEPTH', 'm'),
    MXLR=('NUMBER OF MIXED SIGMA LAYERS', '0--N'),
    OPPT=('OBSERVED PRECIPITATION', 'mm'),
    P10M=('POTENTIAL TEMPERATURE ', 'K'),
    PBLH=('PLANETARY BOUNDARY LAYER HEIGHT', 'm'),
    PRAT=('PRECIPITATION RATE', 'KGM2'),
    PRES=('PRESSURE', 'hPa'),
    PRSS=('SURFACE PRESSURE', 'hPa'),
    PTYP=('PRECIPITATION TYPE (RA=1,TRW=2,ZR=3,ICE=4,SN=5)', '1--5'),
    QSTR=('FLUX MIXING RATIO', '1000.00 ? G/KG'),
    REFC=('COMPOSITE REFLECTIVITY', 'dBZ'),
    RELH=('RELATIVE HUMIDITY', 'PCT'),
    RGHS=('SURFACE ROUGHNESS', 'm'),
    RH2M=('2 METER RELATIVE HUMIDITY', 'PCT'),
    SFCC=('SURFACE CONCENTRATION', '/m3'),
    SHGT=('SURFACE HEIGHT', 'm'),
    SHTF=('SENSIBLE HEAT NET FLUX', 'W/m2'),
    SNOC=('SNOW COVERAGE ', 'PCT'),
    SNOW=('SNOW COVERAGE', '0--1'),
    SNWD=('SNOW DEPTH', 'CM'),
    SOLM=('SOIL MOISTURE', 'kg/m2'),
    SOLT=('SOIL TEMPERATURE ', 'degC'),
    SOLW=('0 TO 200 CM SOIL MOISTURE CONTENT', 'kg/m2'),
    SPHU=('SPECIFIC HUMIDITY', 'kg/kg'),
    STRM=('STREAMLINES ', 'KNTS'),
    T02M=('Temperature at 2 m', 'K'),
    TCLD=('AVERAGE TOTAL CLOUD COVER', 'PCT'),
    TEMP=('Temperature', 'K'),
    THKN=('THICKNESS', '0.10 ? dm'),
    THTS=('SURFACE POTENTIAL TEMPERATURE', 'DEFC'),
    TKEN=('TOTAL KINETIC ENERGY', 'J'),
    TMPS=('Temperature at surface', 'K'),
    TOPO=('TOPOGRAPHY', 'm'),
    TPP1=('1 HOUR ACCUMULATED PRECIPITATION', 'm'),
    TPP3=('3-HOUR ACCUMULATED PRECIPITATION', 'mm'),
    TPP6=('Total precipitation (6-h)', 'm'),
    TPPA=('TOTAL ACCUMULATED PRECIPITATION', 'mm'),
    TPPT=('12-HOUR ACCUMULATED PRECIPITATION', 'mm'),
    TSTR=('FLUX TEMPERATURE', 'degC'),
    U10M=('10 M U-WIND COMPONENT', 'm/s'),
    ULWF=('UPWARD LONG WAVE RADIATION FLUX', 'W/m2'),
    UMOF=('MOMENTUM FLUX, U-WIND COMPONENT', 'N/m2'),
    USTR=('FRICTION VELOCITY', '100.00 ? m/s'),
    UWND=('U-WIND COMPONENT', 'm/s'),
    V10M=('10 M V-WIND COMPONENT', 'm/s'),
    VGTP=('VEGETATION TYPE', '0--1'),
    VMOF=('MOMENTUM FLUX, V-WIND COMPONENT', 'N/m2'),
    VSBY=('VISIBILITY', 'm'),
    VWND=('V-WIND COMPONENT', 'm/s'),
    WESD=('WATER EQUIV. OF ACCUM. SNOW DEPTH', 'kg/m2'),
    WFLX=('WATER VAPOR FLUX', 'kg m2 s'),
    WSPD=('WIND SPEED', 'KNTS'),
    WTMP=('WATER TEMPERATURE', 'degC'),
    WVCT=('WIND VECTORS', 'KNTS'),
    WWND=('W-WIND COMPONENT', 'hPa/s'),
)
# not used in favor of defn on S141.htm
#    HGTS = ('HEIGHT ', '0.10 ? DM'),
#    TMPS = ('SURFACE TEMPERATURE', 'degC'),
#    T02M = ('2 M TEMPERATURE', 'degC'),
#    TEMP = ('TEMPERATURE', 'degC'),
#    TMPS = ('SURFACE TEMPERATURE', 'degC'),
#    TPP6 = ('6-HOUR ACCUMULATED PRECIPITATION', 'mm'),

# not used in favor of defn on https://www.ready.noaa.gov/READYflddesc.php
#    SPHU = ('Specific humidity', 'unknown'),


def readvardef(vheader, out={}):
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
    out['laykeys'] = [(vglvl, keys[vglvl]) for vglvl in vglvls[1:]]
    return out


def writevardef(vglvls, keys, checksums):
    """
    vglvls - iterables of floats [1, .98, ...]
    keys - dictionary of keys by level {vglvli: [vark, ....], ...}
    checksums - integer checsums for each lvl, key --
                {(vglvl, vark): checksum, ...}
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
    fheader = np.fromfile(f, count=1, dtype=thdtype)[0]
    out = {}
    out['LENH'] = hlen = int(fheader['LENH'])
    gridx_off = max(0, (fheader['GRID'][0] - 64) * 1000)
    gridy_off = max(0, (fheader['GRID'][1] - 64) * 1000)
    out['NX'] = int(fheader['NX']) + gridx_off
    out['NY'] = int(fheader['NY']) + gridy_off
    out['NZ'] = int(fheader['NZ'])
    vheader = np.fromfile(f, count=1, dtype='>%dS' % hlen)[0]
    readvardef(vheader, out)
    return out


def getvgtxts(vglvls):
    vgtxts = []
    for vglvl in vglvls:
        if vglvl == 0:
            dp = 0
        else:
            dp = np.floor(np.log10(vglvl) + 1)
        tmpl = '%%6.%df' % np.minimum(5, (5 - dp))
        vgtxt = (tmpl % vglvl)[-6:]
        if len(vgtxt) != 6:
            raise ValueError('Unable to format vglvl:' + str(vglvl))
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

    vglvls = np.append(float(infile.SFCVGLVL),
                       infile.variables['z'][:].array())
    # vgtxts = getvgtxts(vglvls)
    # plus one includes surface
    props['NZ'] = lvar.shape[1] + 1
    props['NY'] = lvar.shape[2]
    props['NX'] = lvar.shape[3]
    datamap = maparlpackedbit(
        path, mode='write', shape=(lvar.shape[0],), props=props)

    theads = datamap['timehead']
    # for key in thdtype.names:
    #    theads[key] = getattr(infile, key)

    # vardefs = datamap['vardef']
    # for ti, vardef in enumerate(vardefs):
    #     sfcvardeftxt = vgtxts[0] + '%2d' % len(svars) +
    #                    ''.join(['%-4s -1 ' % skey.decode()
    #                             for skey in sfckeys])
    #     layvardeftxt = ''
    #     for vgtxt in vgtxts[1:]:
    #         layvardeftxt += vgtxt + '%2d' % len(lvars) +
    #                         ''.join(['%-4s -1 ' % lkey.decode()
    #                                  for lkey in laykeys])
    #
    #
    # vardeftxt = sfcvardeftxt + layvardeftxt
    # defsize = 8+8*len(svars)+(8+8*len(lvars))*len(vgtxts[1:])
    # assert(len(vardeftxt) == defsize)
    # vardefs[:] = vardeftxt
    YYMMDDHHFF = getattr(infile, 'YYMMDDHHFF', '0000000000')
    FF = YYMMDDHHFF[-2:]
    times = gettimes(infile)

    checksums = {}

    for ti, (time, thead) in enumerate(zip(times, theads)):
        for propk in thead.dtype.names:
            if propk in ('NX', 'NY', 'NZ'):
                thead[propk] = '%3d' % props[propk]
            elif propk == 'LENH':
                thead[propk] = '%4d' % datamap['vardef'][ti].itemsize
            else:
                thead[propk] = getattr(infile, propk)
        timestr = time.strftime('%y%m%d%H').encode('ascii') + FF
        thead['YYMMDDHHFF'] = timestr
        for sfck in sfckeys:
            invar = infile.variables[sfck.decode()]
            var_time = datamap['surface'][sfck.decode()]
            varhead = var_time['head']

            _skipprop = ('YYMMDDHHFF', 'LEVEL', 'EXP', 'PREC', 'VAR1')
            for varpropk in varhead.dtype.names:
                if varpropk not in _skipprop:
                    varhead[varpropk][ti] = getattr(invar, varpropk)

            indata = invar[ti]
            CVAR, PREC, NEXP, VAR1, KSUM = pack2d(indata, verbose=False)

            varhead['YYMMDDHHFF'][ti] = timestr
            varhead['LEVEL'][ti] = '%2d' % 0
            varhead['PREC'][ti] = '%14.7E' % PREC
            varhead['EXP'][ti] = '%4d' % NEXP
            varhead['VAR1'][ti] = '%14.7E' % VAR1
            checksums[vglvls[0], sfck] = KSUM
            var_time['data'][ti] = CVAR
        for layk in laykeys:
            invar = infile.variables[layk.decode()]
            var_time = datamap['layers'][layk.decode()][ti]
            for li, var_time_lay in enumerate(var_time):
                varhead = var_time_lay['head']
                for varpropk in varhead.dtype.names:
                    if varpropk not in _skipprop:
                        varhead[varpropk] = getattr(invar, varpropk)

                indata = invar[ti, li]
                CVAR, PREC, NEXP, VAR1, KSUM = pack2d(indata)

                varhead['YYMMDDHHFF'] = timestr
                varhead['LEVEL'] = '%2d' % (li + 1)
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
        thead['LENH'] = datamap['vardef'][ti].itemsize

    datamap.flush()


def maparlpackedbit(path, mode='r', shape=None, props=None):
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
            1 rec of length recl with meta-data
                (thdtype + vardefdtype + hdrdtype)
            for sfckey in sfckeys:
                1 rec of length recl (vhdtype + nx*ny bytes)
            for layer in layers:
                for laykey in laykeys:
                    1 rec of length recl (vhdtype + nx*ny bytes)
    """
    if props is None:
        props = {}

    if props == {}:
        props.update(inqarlpackedbit(path))
    else:
        srflen = 6 + 2 + (4 + 3 + 1) * len(props['sfckeys'])
        laylen = (6 + 2 + (4 + 3 + 1) *
                  len(props['laykeys'])) * (props['NZ'] - 1)
        props['LENH'] = 108 + srflen + laylen

    nx = props['NX']
    ny = props['NY']
    # minus 1 excludes surface
    # nz = props['NZ'] - 1
    hlen = props['LENH']
    sfckeys = props['sfckeys']
    laykeys = props['laykeys']
    ncell = nx * ny
    vardefdtype = dtype('>S%d' % hlen)
    hdrdtype = dtype('>S%d' % (50 + ncell - hlen - thdtype.itemsize))
    lay1dtype = dtype(
        [('head', vhdtype), ('data', dtype('(%d,%d)>1S' % (ny, nx)))])
    sfcdtype = dtype(dict(names=[k.decode() for k in sfckeys], formats=[
                     lay1dtype] * len(sfckeys)))
    layersdtype = dtype([(str(laykey),
                          dtype(dict(names=[k.decode()
                                            for k in layvarkeys],
                                     formats=[lay1dtype] * len(layvarkeys))))
                         for laykey, layvarkeys in laykeys])
    timedtype = dtype([('timehead', thdtype), ('vardef', vardefdtype),
                       ('hdr', hdrdtype), ('surface', sfcdtype),
                       ('layers', layersdtype)])
    datamap = np.memmap(path, timedtype, shape=shape, mode=mode)
    return datamap


def unpack(bytes, VAR1, EXP):
    """
    bytes - nx by ny bytes
    VAR1 - as read directly from LABEL as BYTES with dimension time
    EXP - EXP as read directly from LABEL as BYTES with dimension time
    """
    vold = VAR1.astype('f')
    scale = np.float32(2.0)**np.float32(7 - EXP.astype('i'))
    invscale = np.float32(1.) / scale
    data = (bytes.view('uint8') - np.float32(127.)) * invscale[..., None, None]
    data[..., 0, 0] += vold
    data[..., 0] = np.cumsum(data[..., 0], axis=data.ndim - 2)
    vdata = np.cumsum(data, axis=data.ndim - 1)
    return vdata


def CHAR(c):
    return chr(c).encode('ascii')


def pack2d(RVARA, verbose=False):
    """
    CHARACTER, INTENT(OUT) :: cvar(nxy) ! packed char*1 output array
    REAL,      INTENT(OUT) :: prec      ! precision of packed data array
    INTEGER,   INTENT(OUT) :: nexp      ! packing scaling exponent
    REAL,      INTENT(OUT) :: var1      ! value of real array at position (1,1)
    INTEGER,   INTENT(OUT) :: ksum      ! rotating checksum of packed data
    """
    # MAX = np.maximum
    FLOAT = np.float32
    INT = np.int32
    # ABS = np.abs
    LOG = np.log
    NY, NX = RVARA.shape
    CVAR = np.zeros(RVARA.shape, dtype='uint8')
    RVAR = RVARA.astype('f')
    VAR1 = RVAR[0, 0]

    # find the maximum difference between adjacent elements
    # START ORIGINAL SERIAL CODE
    # ROLD= VAR1
    # RMAX= 0.0
    # for J in range(NY):
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
    colmax = np.abs(np.diff(RVAR, axis=1)).max()
    rowmax = np.abs(np.diff(np.append(VAR1, RVAR[:, 0]), axis=0)).max()
    RMAX = np.maximum(colmax, rowmax)
    # END NUMPY VECTOR CODE

    SEXP = 0.0
    # compute the required scaling exponent
    if RMAX != 0.0:
        SEXP = LOG(RMAX) / LOG(np.float32(2.))

    NEXP = INT(SEXP)
    # positive or whole number scaling round up for lower precision
    if SEXP >= 0.0 or (SEXP % 1.0) == 0.0:
        NEXP = NEXP + 1
    # precision range is -127 to 127 or 254
    PREC = np.float32((2.0**NEXP) / 254.0)
    SCEXP = np.float32(2.0**(7 - NEXP))

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
        ROLD = FLOAT(ICVAL - 127) / SCEXP + ROLD
        ROLDS[J] = ROLD

    ROLD = ROLDS
    for I in range(1, NX):
        ICVAL = INT((RVAR[:, I] - ROLD) * SCEXP + 127.5)
        CVAR[:, I] = ICVAL
        ROLD = FLOAT(ICVAL - 127) / SCEXP + ROLD
    KSUM = INT(CVAR.sum()) % 255
    # END NUMPY VECTOR CODE

    # assert((CVART == CVAR).all())
    # assert(KSUM == KSUMT)
    return CVAR.view('>S1'), PREC, NEXP, VAR1, KSUM


class arlpackedbit(PseudoNetCDFFile):
    """
    arlpackedbit reads files formatted according to NOAA's arl packed bit
    format

    Format as follows:

    for t in times:
        thdtype = dtype([('YYMMDDHHFF', '>10S'), ('LEVEL', '>2S'),
                         ('GRID', '>2S') , ('INDX', '>4S'), ('Z1', '>4S'),
                         ('MB', '2>14S'),  ('SRCE', '>S4'), ('MGRID', '>S3'),
                         ('VSYS', '>S2'), ('POLLAT', '>S7'), ('POLLON', '>S7'),
                         ('REFLAT','>S7'), ('REFLON','>S7'), ('GRIDX','>S7'),
                         ('ORIENT','>S7'), ('TANLAT','>S7'), ('SYNCHX','>S7'),
                         ('SYNCHY','>S7'), ('SYNCHLAT','>S7'),
                         ('SYNCHLON','>S7'), ('RESERVED', '>S7'),
                         ('NX', '>S3'), ('NY', '>S3'), ('NZ', '>S3'),
                         ('VSYS2', '>S2'), ('LENH', '>S4')]) +
                         ny * nx - 158 - 50
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
    P(ny,nx) = ny by nx 1-byte values dervied from real values (R) according
               to the formula below.

    P(i,j) = (Ri,j  - Ri-1,j)* (2**(7-(ln dRmax / ln 2)))
    """
    @classmethod
    def isMine(cls, path):
        testchunk = np.fromfile(path,
                                dtype=dtype([('YYMMDDHHFF', '>10S'),
                                             ('LEVEL', '>2S'), ('GRID', '>2S'),
                                             ('INDX', '>4S')]), count=1)
        check = testchunk[0]['INDX'] == b'INDX'
        return check

    def __init__(
        self, path, shape=None, cache=False,
        earth_radius=6370000, synchlon=None, synchlat=None
    ):
        """
        Parameters
        ----------
        path : string
            path to arl packed bit formatted file
        shape : tuple or None
            shape of file to over ride derived shape
        earth_radius : scalar
            radius of the earth used in converting projected to lat/lon
        synchlon: scalar
            The SYNCHLON variable has 6 digits of precision in the file
            this keyword allows you to provide more.
        synchlat: scalar
            The SYNCHLAT variable has 6 digits of precision in the file
            this keyword allows you to provide more.

        Returns
        -------
        arlf : arlpackedbit
            PseudoNetCDF file with packed bit contents unpacked
        """
        self._path = path
        self._f = f = open(path, 'r')
        self._cache = cache
        self._earth_radius = earth_radius
        # f.seek(0, 2)
        # fsize = f.tell()
        f.seek(0, 0)

        # vord = np.vectorize(ord)
#        timesfc = []
#        times = []
#        tflag = []

        f.seek(0, 0)
        props = {}
        self._datamap = maparlpackedbit(path, shape=shape, props=props)
        datamap = self._datamap
        t0hdr = datamap['timehead'][0]
        for key in thdtype.names:
            setattr(self, key, t0hdr[key])

        out = readvardef(datamap['vardef'][0])
        alllayvarkeys = []
        for layk, layvarkeys in out['laykeys']:
            alllayvarkeys.extend(
                [k.decode() for k in layvarkeys
                 if k.decode() not in alllayvarkeys])
        self._layvarkeys = tuple(alllayvarkeys)
        tflag = np.char.replace(datamap['timehead']['YYMMDDHHFF'], b' ', b'0')
        sfckeys = self._datamap['surface'].dtype.names
        laykeys = self._layvarkeys

        self.variables = PseudoNetCDFVariables(
            func=self._getvar, keys=list(sfckeys + laykeys))
        self.createDimension('time', datamap.shape[0])
        # minus 1 excludes surface
        self.createDimension('z', int(self.NZ) - 1)
        if synchlon is None:
            self._synchlon = float(self.SYNCHLON)
        else:
            self._synchlon = synchlon

        if synchlat is None:
            self._synchlat = float(self.SYNCHLAT)
        else:
            self._synchlat = synchlat

        gridx = float(t0hdr['GRIDX']) * 1000.
        nx = props['NX']
        ny = props['NY']
        if gridx != 0:
            self.XSIZE = self.YSIZE = gridx
            self.createDimension('x', nx)
            self.createDimension('y', ny)
            x = self.createVariable('x', 'f', ('x',))
            x[:] = np.arange(x[:].size) * gridx
            y = self.createVariable('y', 'f', ('y',))
            y[:] = np.arange(y[:].size) * gridx
            self._addcrs()
        else:
            self.createDimension('x', nx)
            self.createDimension('y', ny)
            self.createDimension('nv', 2)
            x = self.createVariable('x', 'f', ('x',))
            x.standard_name = 'longitude'
            x.units = 'degrees'
            xe = self.createVariable('x_bounds', 'f', ('x', 'nv'))
            xe.standard_name = 'longitude_bounds'
            xe.units = 'degrees'
            x[:] = np.arange(0, nx) * float(self.REFLON) + self._synchlon
            xe[:-1, 1] = x[:1] + np.diff(x) / 2
            xe[1:, 0] = x[1:] - np.diff(x) / 2
            xe[0, 0] = x[0] - np.diff(x)[0] / 2
            xe[1, 1] = x[1] + np.diff(x)[1] / 2
            y = self.createVariable('y', 'f', ('y',))
            y.standard_name = 'latitude'
            y.units = 'degrees'
            ye = self.createVariable('y_bounds', 'f', ('y', 'nv'))
            ye.standard_name = 'latitude_bounds'
            ye.units = 'degrees'
            y[:] = np.arange(0, ny) * float(self.REFLAT) + self._synchlat
            ye[:-1, 1] = y[:1] + np.diff(y) / 2
            ye[1:, 0] = y[1:] - np.diff(y) / 2
            ye[0, 0] = y[0] - np.diff(y)[0] / 2
            ye[1, 1] = y[1] + np.diff(y)[1] / 2

        times = [datetime.strptime(
            t.astype('S8').decode(), '%y%m%d%H') for t in tflag]
        rtime = times[0]
        hours_since = [(t - rtime).total_seconds() // 3600 for t in times]
        timev = self.createVariable('time', 'i', ('time',))
        timev[:] = hours_since
        timev.units = rtime.strftime('hours since %F %H:%M:%S')
        z = self.createVariable('z', 'f', ('z',))
        z[:] = out['vglvls'][1:]
        self.SFCVGLVL = out['vglvls'][0]
        _units = {1: 'pressure sigma', 2: 'pressure absolute',
                  3: 'terrain sigma', 4: 'hybrid sigma'}
        z.units = _units.get(int(self.VSYS2), 'unknown')
        self.setCoords(['time', 'z', 'y', 'x'])

    def _addcrs(self):
        x = self.variables['x']
        y = self.variables['y']
        gridx = np.diff(x[:])[0]
        gridy = np.diff(y[:])[0]
        nx = x.size
        ny = y.size
        crs = self.createVariable('crs', 'i', ())
        s = self
        tanlat = np.float64(s.TANLAT)
        reflon = np.float64(s.REFLON)
        reflat = np.float64(s.REFLAT)
        atanlat = np.abs(tanlat)
        pname = crs.grid_mapping_name = {
            0: 'equatorial_mercator',
            90: 'polar_stereographic',
        }.get(atanlat, 'lambert_conformal_conic')
        if pname == 'lambert_conformal_conic':
            crs.standard_parallel = np.array([tanlat, tanlat])
            crs.longitude_of_central_meridian = reflon
            crs.latitude_of_projection_origin = reflat
        elif pname == 'polar_stereographic':
            crs.straight_vertical_longitude_from_pole = reflon
            crs.standard_parallel = reflat
            crs.latitude_of_projection_origin = tanlat
        else:
            raise KeyError('Not yet implemented equatorial mercator')

        crs.earth_radius = self._earth_radius  # WRF-based radius
        scellx = (float(self.SYNCHX) - 1) * gridx  # 1-based I cell
        scelly = (float(self.SYNCHY) - 1) * gridy  # 1-based J cell
        slon, slat = self._synchlon, self._synchlat
        if slon > 180:
            slon = slon % 180 - 180

        halfwidth = gridx * (nx - 1) / 2
        halfheight = gridy * (ny - 1) / 2
        crs.false_easting = halfwidth
        crs.false_northing = halfheight
        llcrnrlon, llcrnrlat = self.xy2ll(scellx, scelly)

        x_within_precision = np.round(slon / llcrnrlon, 3) == 1
        y_within_precision = np.round(slat / llcrnrlat, 4) == 1
        if not (x_within_precision and y_within_precision):
            warn(
                'Grid not centered; using SYNCHLAT/SYNCHLON ' +
                'with limited to 6 significant digits ' +
                'to calculate false easting/northing'
            )
            crs.false_easting = 0.
            crs.false_northing = 0.
            llcrnrx, llcrnry = self.ll2xy(slon, slat)
            crs.false_easting = -llcrnrx + scellx
            crs.false_northing = -llcrnry + scelly

        self.Conventions = 'CF-1.6'

    def _getvar(self, varkey):
        datamap = self._datamap
        stdname, stdunit = stdprops.get(varkey, (varkey, 'unknown'))
        if varkey in datamap['surface'].dtype.names:
            vhead = datamap['surface'][varkey]['head']
            v11 = vhead['VAR1']
            EXP = vhead['EXP']
            props = dict([(pk, vhead[pk][0]) for pk in vhead.dtype.names
                          if pk not in ('YYMMDDHHFF', 'LEVEL')])
            bytes = datamap['surface'][varkey]['data']
            vdata = unpack(bytes, v11, EXP)
            # if varkey == 'MXLR':
            #  CVAR, PREC, NEXP, VAR1, KSUM = pack2d(vdata[21], verbose = True)

            out = PseudoNetCDFVariable(
                self, varkey, 'f', ('time', 'y', 'x'),
                values=vdata, units=stdunit, standard_name=stdname, **props
            )
        elif varkey in self._layvarkeys:
            laykeys = datamap['layers'].dtype.names
            mylaykeys = [laykey for laykey in laykeys
                         if varkey in datamap['layers'][laykey].dtype.names]
            bytes = np.array([datamap['layers'][lk][varkey]['data']
                              for lk in mylaykeys]).swapaxes(0, 1)
            vhead = np.array([datamap['layers'][lk][varkey]['head']
                              for lk in mylaykeys]).swapaxes(0, 1)
            v11 = vhead['VAR1']
            EXP = vhead['EXP']
            props = dict([(pk, vhead[pk][0, 0]) for pk in vhead.dtype.names
                          if pk not in ('YYMMDDHHFF', 'LEVEL')])
            props['LEVEL_START'] = vhead['LEVEL'][0, 0]
            props['LEVEL_END'] = vhead['LEVEL'][-1, -1]
            vdata = unpack(bytes, v11, EXP)
            out = PseudoNetCDFVariable(
                self, varkey, 'f', ('time', 'z', 'y', 'x'),
                values=vdata, units=stdunit,
                standard_name=stdname, **props
            )
        if self._cache:
            self.variables[varkey] = out
        return out

    def getMap(self, maptype='basemap_auto', **kwds):
        if 'latitude_bounds' in self.variables:
            return PseudoNetCDFFile.getMap(self, maptype=maptype, **kwds)

        from PseudoNetCDF.coordutil import basemap_from_proj4
        kwds = kwds.copy()
        myproj = self.getproj(withgrid=True, projformat='pyproj')
        myprojstr = self.getproj(withgrid=True, projformat='proj4')
        llcrnrlon, llcrnrlat = myproj(
            self.variables['x'][0], self.variables['y'][0], inverse=True
        )
        urcrnrlon, urcrnrlat = myproj(
            self.variables['x'][-1], self.variables['y'][-1], inverse=True
        )
        kwds['llcrnrlon'] = llcrnrlon
        kwds['llcrnrlat'] = llcrnrlat
        kwds['urcrnrlon'] = urcrnrlon
        kwds['urcrnrlat'] = urcrnrlat
        return basemap_from_proj4(myprojstr, **kwds)

    def plot(self, *args, **kwds):
        kwds.setdefault('plottype', 'x-y')
        ax = PseudoNetCDFFile.plot(self, *args, **kwds)
        if kwds['plottype'] == 'x-y':
            map_kw = kwds.get('map_kw', None)
            if map_kw is None:
                map_kw = {}
            try:
                map_kw = map_kw.copy()
                coastlines = map_kw.pop('coastlines', True)
                countries = map_kw.pop('countries', True)
                states = map_kw.pop('states', False)
                counties = map_kw.pop('counties', False)
                bmap = self.getMap(**map_kw)
                if coastlines:
                    bmap.drawcoastlines(ax=ax)
                if countries:
                    bmap.drawcountries(ax=ax)
                if states:
                    bmap.drawstates(ax=ax)
                if counties:
                    bmap.drawcounties(ax=ax)
            except Exception:
                pass
        return ax

    def getproj(self, withgrid=False, projformat='pyproj'):
        from PseudoNetCDF.coordutil import getproj4_from_cf_var
        gridmapping = self.variables['crs']
        proj4str = getproj4_from_cf_var(gridmapping, withgrid=withgrid)
        preserve_units = withgrid
        if projformat == 'proj4':
            return proj4str
        elif projformat == 'pyproj':
            import pyproj
            # pyproj adds +units=m, which is not right for latlon/lonlat
            if '+proj=lonlat' in proj4str or '+proj=latlon' in proj4str:
                preserve_units = True
            return pyproj.Proj(proj4str, preserve_units=preserve_units)
        elif projformat == 'wkt':
            import osr
            srs = osr.SpatialReference()
            # Imports WKT to Spatial Reference Object
            srs.ImportFromProj4(proj4str)
            srs.ExportToWkt()  # converts the WKT to an ESRI-compatible format
            return srs.ExportToWkt()
        else:
            raise ValueError('projformat must be pyproj, proj4 or wkt')


registerwriter('noaafiles.arlpackedbit', writearlpackedbit)
registerwriter('arlpackedbit', writearlpackedbit)

if __name__ == '__main__':
    import sys
    out = arlpackedbit(sys.argv[1])
