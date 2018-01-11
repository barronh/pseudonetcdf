#!/usr/bin/env python
from __future__ import print_function

import sys
import os
import re
import warnings

from collections import defaultdict
import textwrap
import json
from datetime import datetime, timedelta


import numpy as np
from numpy import exp, log, pi, zeros, nan

from netCDF4 import Dataset
netcdf = Dataset
from PseudoNetCDF import getvarpnc, slice_dim, extract
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF.pncgen import pncgen
from PseudoNetCDF.coordutil import gettimes
from PseudoNetCDF.geoschemfiles import bpch
from PseudoNetCDF.register import registerreader

messages = ""
def formatwarning(message, category, filename, lineno, line = 0):
    """
    Warning formatter function
    """
    global messages
    strout = "\n***********\n%s:%s: %s:\n\t%s\n***********\n" % (filename, lineno, category.__name__, message)
    messages += strout
    return strout

warnings.formatwarning = formatwarning

warn = warnings.warn

def interpbound(x, xp, fp):
    "Standard numpy interp with minimum=fp.min and max=fp.max"
    sfc = fp[0]
    top = fp[-1]
    out = np.interp(x[::-1], xp[::-1], fp[::-1], left = top, right = sfc)[::-1]
    return out


_GCLIST = "NO2 NO O3 NO3 OH HO2 N2O5 HNO3 HONO PNA H2O2 NTR ROOH FORM ALD2 PAR CO MEPX FACD C2O3 PAN PACD AACD PANX OLE ETH IOLE TOL CRES OPEN MGLY XYL ISOP SO2 SULF ETHA BENZENE NH3 SV_ALK SV_XYL1 SV_XYL2 SV_TOL1 SV_TOL2 SV_BNZ1 SV_BNZ2 SV_TRP1 SV_TRP2 SV_ISO1 SV_ISO2 SV_SQT HG HGIIGAS MACR MVK".split()
_NRLIST = []
_AELIST = "ASO4J ASO4I AALKJ AXYL1J AXYL2J AXYL3J ATOL1J ATOL2J ATOL3J ABNZ1J ABNZ2J ABNZ3J ATRP1J ATRP2J AISO1J AISO2J ASQTJ ACORS ASOIL AISO3J AOLGAJ          AOLGBJ ANIJ ACR_IIIJ ACR_VIJ APBJ APBK ACDJ AMN_HAPSJ AMN_HAPSK APHGJ AORGPAJ AORGPAI".split()
_NMLIST  = "NUMATKN NUMACC NUMCOR".split()
_SFLIST = "SRFATKN SRFACC SRFCOR".split()

_CBSPCS = r"""
        "O3": {
            "expression": "O3",
            "outunit": "ppmV"
        },
        "NO": {
            "expression": "NO",
            "outunit": "ppmV"
        },
        "NO2": {
            "expression": "NO2",
            "outunit": "ppmV"
        },
        "NO3": {
            "expression": "NO3",
            "outunit": "ppmV"
        },
        "NH3": {
            "expression": "NH3",
            "outunit": "ppmV"
        },
        "N2O5": {
            "expression": "N2O5",
            "outunit": "ppmV"
        },
        "HONO": {
            "expression": "HNO2",
            "outunit": "ppmV"
        },
        "PNA": {
            "expression": "HNO4",
            "outunit": "ppmV"
        },
        "SO2": {
            "expression": "SO2",
            "outunit": "ppmV"
        },
        "CO": {
            "expression": "CO",
            "outunit": "ppmV"
        },
        "FORM": {
            "expression": "CH2O",
            "outunit": "ppmV"
        },
        "ALD2": {
            "expression": "ALD2",
            "outunit": "ppmV"
        },
        "ALDX": {
            "expression": "RCHO",
            "outunit": "ppmV"
        },
        "OLE": {
            "expression": "0.5 * PRPE",
            "outunit": "ppmV",
            "comment1": "Half the carbon is assumed to be propene and half butene",
            "comment2": "propene maps to 1 PAR and 1 OLE."
        },
        "IOLE": {
            "expression": "0.5 * 1./4. * 3. * PRPE",
            "outunit": "ppmV",
            "comment1": "Half the carbon is assumed to be propene and half butene",
            "comment2": "IOLE is a 4 carbon species.",
            "comment2": "PRPE is a 3 carbon species."
        },
        "ETHA": {
            "expression": "C2H6",
            "outunit": "ppmV"
        },
        "TOL": {
            "expression": "TOLU",
            "outunit": "ppmV"
        },
        "XYL": {
            "expression": "XYLE",
            "outunit": "ppmV"
        },
        "ISOP": {
            "expression": "ISOP",
            "outunit": "ppmV"
        },
        "ISPD": {
            "expression": "MACR + MVK",
            "outunit": "ppmV"
        },
        "HNO3": {
            "expression": "HNO3",
            "outunit": "ppmV"
        },
        "NTR": {
            "expression": "R4N2",
            "outunit": "ppmV"
        },
        "H2O2": {
            "expression": "H2O2",
            "outunit": "ppmV"
        },
        "MEPX": {
            "expression": "MP",
            "outunit": "ppmV"
        },
        "PAN": {
            "expression": "PAN",
            "outunit": "ppmV"
        },
        "PANX": {
            "expression": "PPN + PMN",
            "outunit": "ppmV"
        },
        "PACD": {
            "expression": "MAP",
            "outunit": "ppmV"
        },"""

_CB05SPCS = r"""
        "PAR": {
            "expression": "1.5 * C3H8 + 4. * ALK4 + 3. * ACET + 4. * MEK + 1. * BENZ + 0.5 * PRPE",
            "outunit": "ppmV",
            "comment1": "Propane maps as 1.5 PAR and 1.5 unreactive",
            "comment2": "ALK4 is a 4C species in GEOS-Chem and maps as 4 PAR",
            "comment3": "Official mapping of acetone for CB05 is 3 PAR, but should likely be 2 or less",
            "comment4": "Methyl ethyl ketone maps a 4 PAR.",
            "comment5": "Benzene maps as 1PAR and 1BENZENE for SOA",
            "comment6": "Half the carbon is assumed to be propene which is 1 PAR + 1 OLE"
        },
        "ROOH": {
            "expression": "RIP",
            "outunit": "ppmV"
        },"""

_CB6SPCS = r"""
        "PAR": {
            "expression": "4. * ALK4 + 0.5 * PRPE",
            "outunit": "ppmV",
            "comment1": "Propane maps as 1.5 PAR and 1.5 unreactive",
            "comment2": "ALK4 is a 4C species in GEOS-Chem and maps as 4 PAR",
            "comment5": "Missing; Benzene maps as 1PAR and 1BENZENE for SOA",
            "comment6": "Half the carbon is assumed to be propene which is 1 PAR + 1 OLE"
        },
        "GLYD": {
            "expression": "GLYC",
            "outunit": "ppmV",
            "comment1": "glycoaldehyde from both models"
        },
        "EPOX": {
            "expression": "IEPOX",
            "outunit": "ppmV",
            "comment1": "EPOX and IEPOX are both isoprene epoxides"
        },
        "INTR": {"expression": "ISOPN",
            "outunit": "ppmV",
            "comment1": "INTR is the isoprene nitrate and so is ISOPN"
        },
        "ISPX": {
            "expression": "RIP",
            "outunit": "ppmV",
            "comment1": "RIP is isoprene based higher peroxide"
        },"""

_GC_TO_AE6_COMMON = r"""
        "AALJ": {
            "expression": "0.05695 * DST1",
            "outunit": "micrograms/m**3"
        },
        "AECI": {
            "expression": "0.1 * BCPI + 0.1 * BCPO",
            "outunit": "micrograms/m**3"
        },
        "AECJ": {
            "expression": "0.9 * BCPI + 0.9 * BCPO",
            "outunit": "micrograms/m**3"
        },
        "APOCI": {
            "expression": "0.1 * OCPI + 0.1 * OCPO",
            "outunit": "micrograms/m**3"
        },
        "APOCJ": {
            "expression": "0.9 * OCPI + 0.9 * OCPO + 0.01075 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ACAJ": {
            "expression": "0.0118 * SALA + 0.07940 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ACLJ": {
            "expression": "0.00945 * DST1 + 0.5538 * SALA",
            "outunit": "micrograms/m**3"
        },
        "ACLK": {
            "expression": "0.01190 * DST2 + 0.01190 * DST3 + 0.01190 * DST4 + 0.5538 * SALC",
            "outunit": "micrograms/m**3"
        },
        "AFEJ": {
            "expression": "0.03355 * DST1",
            "outunit": "micrograms/m**3"
        },
        "AKJ": {
            "expression": "0.0114 * SALA + 0.03770 * DST1",
            "outunit": "micrograms/m**3"
        },
        "AMGJ": {
            "expression": "0.0368 * SALA",
            "outunit": "micrograms/m**3"
        },
        "AMNJ": {
            "expression": "0.00115 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ANAJ": {
            "expression": "0.3086 * SALA + 0.03935 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ANH4I": {
            "expression": "0.01 * NH4",
            "outunit": "micrograms/m**3"
        },
        "ANH4J": {
            "expression": "0.00005 * DST1 + 0.99 * NH4",
            "outunit": "micrograms/m**3"
        },
        "ANO3I": {
            "expression": "0.01 * NIT",
            "outunit": "micrograms/m**3"
        },
        "ANO3J": {
            "expression": "0.00020 * DST1 + 0.99 * NIT",
            "outunit": "micrograms/m**3"
        },
        "ANO3K": {
            "expression": "0.0016 * DST2 + 0.0016 * DST3 + 0.0016 * DST4 + NITs",
            "outunit": "micrograms/m**3"
        },
        "AOTHRJ": {
            "expression": "0.50219 * DST1",
            "outunit": "micrograms/m**3"
        },
        "APNCOMI": {
            "expression": "0.4 * 0.1 * OCPI + 0.4 * 0.1 * OCPO",
            "outunit": "micrograms/m**3"
        },
        "APNCOMJ": {
            "expression": "0.0043 * DST1 + 0.4 * 0.9 * OCPI + 0.4 * 0.9 * OCPO",
            "outunit": "micrograms/m**3"
        },
        "ASEACAT": {
            "expression": "0.3685 * SALC",
            "outunit": "micrograms/m**3"
        },
        "ASIJ": {
            "expression": "0.19435 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ASO4I": {
            "expression": "0.01 * SO4",
            "outunit": "micrograms/m**3"
        },
        "ASO4J": {
            "expression": "0.99 * SO4 + 0.0225 * DST1 + 0.0776 * SALA",
            "outunit": "micrograms/m**3"
        },
        "ASO4K": {
            "expression": "0.0776 * SALC + 0.02655 * DST2 + 0.02655 * DST3 + 0.02655 * DST4 + SO4s",
            "outunit": "micrograms/m**3"
        },
        "ASOIL": {
            "expression": "0.95995 * DST2 + 0.95995 * DST3 + 0.95995 * DST4",
            "outunit": "micrograms/m**3"
        },
        "ATIJ": {
            "expression": "0.0028 * DST1",
            "outunit": "micrograms/m**3"
        },"""

_GC8_to_AE6 = r"""
        "ABNZ1J": {
            "expression": "0.12 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "ABNZ2J": {
            "expression": "0.04 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "ABNZ3J": {
            "expression": "0.32 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "AISO1J": {
            "expression": "0.75 * SOA4",
            "outunit": "micrograms/m**3"
        },
        "AISO2J": {
            "expression": "0.25 * SOA4",
            "outunit": "micrograms/m**3"
        },
        "ASQTJ": {
            "expression": "SOA3",
            "outunit": "micrograms/m**3"
        },
        "ATOL1J": {
            "expression": "0.04 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "ATOL2J": {
            "expression": "0.04 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "ATOL3J": {
            "expression": "0.29 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "ATRP1J": {
            "expression": "0.33 * SOA1",
            "outunit": "micrograms/m**3"
        },
        "ATRP1J": {
            "expression": "0.33 * SOA2",
            "outunit": "micrograms/m**3"
        },
        "ATRP2J": {
            "expression": "0.67 * SOA1",
            "outunit": "micrograms/m**3"
        },
        "ATRP2J": {
            "expression": "0.67 * SOA2",
            "outunit": "micrograms/m**3"
        },
        "AXYL1J": {
            "expression": "0.03 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "AXYL2J": {
            "expression": "0.01 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "AXYL3J": {
            "expression": "0.11 * SOA5",
            "outunit": "micrograms/m**3"
        },
        "SV_BNZ1": {
            "expression": "0.06 * SOG5",
            "outunit": "ppmV"
        },
        "SV_BNZ2": {
            "expression": "0.23 * SOG5",
            "outunit": "ppmV"
        },
        "SV_ISO1": {
            "expression": "0.75 * SOG4",
            "outunit": "ppmV"
        },
        "SV_ISO2": {
            "expression": "0.25 * SOG4",
            "outunit": "ppmV"
        },
        "SV_SQT": {
            "expression": "SOG3",
            "outunit": "ppmV"
        },
        "SV_TOL1": {
            "expression": "0.23 * SOG5",
            "outunit": "ppmV"
        },
        "SV_TOL2": {
            "expression": "0.23 * SOG5",
            "outunit": "ppmV"
        },
        "SV_TRP1": {
            "expression": "0.33 * SOG1",
            "outunit": "ppmV"
        },
        "SV_TRP1": {
            "expression": "0.33 * SOG2",
            "outunit": "ppmV"
        },
        "SV_TRP2": {
            "expression": "0.67 * SOG1",
            "outunit": "ppmV"
        },
        "SV_TRP2": {
            "expression": "0.67 * SOG2",
            "outunit": "ppmV"
        },
        "SV_XYL1": {
            "expression": "0.19 * SOG5",
            "outunit": "ppmV"
        },
        "SV_XYL2": {
            "expression": "0.06 * SOG5",
            "outunit": "ppmV"
        }"""

_GCgt8_to_AE6 = r"""
        "BENZENE": {
            "expression": "BENZ",
            "outunit": "ppmV"
        },
        "ABNZ1J": {
            "expression": "0.12 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ABNZ2J": {
            "expression": "0.04 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ABNZ3J": {
            "expression": "0.32 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "AISO1J": {
            "expression": "0.75 * (ISOA1 + ISOA2 + ISOA3)",
            "outunit": "micrograms/m**3"
        },
        "AISO2J": {
            "expression": "0.25 * (ISOA1 + ISOA2 + ISOA3)",
            "outunit": "micrograms/m**3"
        },
        "ASQTJ": {
            "expression": "0.33 * TSOA0 + 0.33 * TSOA1 + 0.33 * TSOA2 + 0.33 * TSOA3",
            "outunit": "micrograms/m**3"
        },
        "ATOL1J": {
            "expression": "0.04 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ATOL2J": {
            "expression": "0.04 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ATOL3J": {
            "expression": "0.29 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ATRP1J": {
            "expression": "0.33 * TSOA0 + 0.33 * TSOA1 + 0.33 * TSOA2 + 0.33 * TSOA3",
            "outunit": "micrograms/m**3"
        },
        "ATRP2J": {
            "expression": "0.34 * TSOA0 + 0.34 * TSOA1 + 0.34 * TSOA2 + 0.34 * TSOA3",
            "outunit": "micrograms/m**3"
        },
        "AXYL1J": {
            "expression": "0.03 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "AXYL2J": {
            "expression": "0.01 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "AXYL3J": {
            "expression": "0.11 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "SV_BNZ1": {
            "expression": "0.06 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_BNZ2": {
            "expression": "0.23 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_ISO1": {
            "expression": "0.75 * (ISOG1 + ISOG2 + ISOG3)",
            "outunit": "ppmV"
        },
        "SV_ISO2": {
            "expression": "0.25 * (ISOG1 + ISOG2 + ISOG3)",
            "outunit": "ppmV"
        },
        "SV_SQT": {
            "expression": "0.33 * TSOG0 + 0.33 * TSOG1 + 0.33 * TSOG2 + 0.33 * TSOG3",
            "outunit": "ppmV"
        },
        "SV_TOL1": {
            "expression": "0.23 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_TOL2": {
            "expression": "0.23 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_TRP1": {
            "expression": "0.33 * TSOG0 + 0.33 * TSOG1 + 0.33 * TSOG2 + 0.33 * TSOG3",
            "outunit": "ppmV"
        },
        "SV_TRP2": {
            "expression": "0.34 * TSOG0 + 0.34 * TSOG1 + 0.34 * TSOG2 + 0.34 * TSOG3",
            "outunit": "ppmV"
        },
        "SV_XYL1": {
            "expression": "0.19 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_XYL2": {
            "expression": "0.06 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        }"""

def geticbcnames(inpath):
    import io
    import re
    import numpy as np
    with open(inpath) as infile:
        intxt = infile.read()
    result = re.findall("(TYPE_\S+\s*=\s*(('[^\n]+',\s*)+))", intxt)
    names = result[0][1].strip()[1:-1].replace("'", '') + '\n'
    nflds = len(names.split(':'))
    if nflds == 16: fmts = 'S16 f c f c f c f S16 f c c b b b b'.split()
    else: fmts = 'S16 f c f c f c f S16 f c b b b b'.split()
    data = re.sub('\s+', '', result[1][1])[:-1]
    data = re.sub("'", '', data)
    data = re.sub(",", '\n', data)
    data = (names + data).encode()
    rdata = np.recfromcsv(io.BytesIO(data), delimiter = b':', dtype = fmts)
    icbc_sur = []
    for row in rdata:
        if row['icbc_sur'] != b'': icbc_sur.append(row['icbc_sur'].decode().strip())
    icbc_sur =  set(icbc_sur)
    result = icbc_sur.union([s.decode().strip() for s in rdata['spc']])
    return list(result)

def getspclists(args):
    if args.GCNML is None:
        GCLIST = _GCLIST
    else:
        GCLIST = geticbcnames(args.GCNML)
    if args.AENML is None:
        AELIST = _AELIST
        NMLIST  = _NMLIST
        SFLIST = _SFLIST
    else:
        AELIST = geticbcnames(args.AENML)
        SFLIST = [s for s in AELIST if s.startswith('SRF')]
        NMLIST = [s for s in AELIST if s.startswith('NUM')]
        AELIST = [s for s in AELIST if s not in (SFLIST + NMLIST)]
    if args.NRNML is None:
        NRLIST = _NRLIST
    else:
        NRLIST = geticbcnames(args.NRNML)
    GCLIST += NRLIST
    return dict(GCLIST = GCLIST, AELIST = AELIST, SFLIST = SFLIST, NMLIST = NMLIST)

def makedefaulticon(metcroprops, args):
    """
    Create an empty bcon file whose dimensions
    are consistent with metcroprops
    
    metcroprops - dictionary of properties from IOAPI metadata
    """
    profilepath = args.ICONPROFILEPATH
    SPCLISTS = getspclists(args)
    GCLIST = SPCLISTS['GCLIST']
    AELIST = SPCLISTS['AELIST']
    SFLIST = SPCLISTS['SFLIST']
    NMLIST = SPCLISTS['NMLIST']
    if not profilepath is None:
        from PseudoNetCDF.cmaqfiles import icon_profile
        iconfile = icon_profile(profilepath)
        if not args.keepall: iconfile = getvarpnc(iconfile, GCLIST + AELIST + SFLIST + NMLIST)
        varlist = [varkey for varkey in iconfile.variables.keys() if varkey in (GCLIST + AELIST + SFLIST + NMLIST)]
        missingunits = [varkey for varkey, var in iconfile.variables.items() if var[:].ndim == 2 and varkey not in GCLIST + AELIST + NMLIST + SFLIST]
        missing = ', '.join(missingunits)
        GCLIST = missingunits + GCLIST
        if missing != '': warn(profilepath + ' contains ' + missing + ' whose units are assumed ppmV')
        def getvals(varkey):
            if varkey not in iconfile.variables:
                warn(varkey + ' not found in icon profile path')
                return 1e-32
            invglvls = np.convolve([.5, .5], metcroprops['VGLVLS'], mode = 'valid')
            outvglvls = np.convolve([.5, .5], iconfile.VGLVLS, mode = 'valid')
            invar = iconfile.variables[varkey]
            allvals = interpbound(invglvls, outvglvls, invar[:])[:, None, None].repeat(metcroprops['NROWS'], 1).repeat(metcroprops['NCOLS'], 2)
            return allvals
    else:
        def getvals(varkey):
            return 1e-32
 
    DI = Dataset('dummyicon.nc', 'w', format = 'NETCDF3_CLASSIC')
    DI.createDimension('TSTEP', None)
    DI.createDimension('DATE-TIME', 2) ;
    DI.createDimension('VAR', 106) ;
    DI.createDimension('LAY', metcroprops['NLAYS']) ;
    DI.createDimension('ROW', metcroprops['NROWS']);
    DI.createDimension('COL', metcroprops['NCOLS']) ;

    TFLAG = DI.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME')) ;
    TFLAG.units = "<YYYYDDD,HHMMSS>" ;
    TFLAG.long_name = "TFLAG           " ;
    TFLAG.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
    AIRDEN = DI.createVariable('AIRDEN', 'f', ('TSTEP', 'LAY', 'ROW', 'COL')) ;
    AIRDEN.long_name = "AIRDEN          " ;
    AIRDEN.units = "molec/cm3       " ;
    AIRDEN.coordinates = "lon lat" ;
    AIRDEN.var_desc = "AIRDEN          " ;
    def getunit(varkey):
        if varkey in GCLIST:
            return 'ppmV'
        elif varkey in AELIST:
            return 'micrograms/m**3'
        elif varkey in NMLIST:
            return '#/m**3'
        elif varkey in SFLIST:
            return 'm**2/m**3'
    
    #varkey_unit = [(i, 'ppmV') for i in GCLIST] + [(i, 'micrograms/m**3') for i in AELIST] + [(i, '#/m**3') for i in NMLIST] + [(i, 'm**2/m**3') for i in SFLIST]
    varkey_unit = [(varkey, getunit(varkey)) for varkey in varlist]
    for varkey, unit in varkey_unit:
        var = DI.createVariable(varkey, 'f', ('TSTEP', 'LAY', 'ROW', 'COL'))
        var.long_name = varkey.ljust(16);
        var.units = unit.ljust(16);
        var.coordinates = "lon lat" ;
        var.var_desc = ("Variable %s" % varkey).ljust(80);
    
    for varkey, unit in varkey_unit:
        DI.variables[varkey][0] = getvals(varkey)
    
    TFLAG[0] = 0
    DI.IOAPI_VERSION = "$Id: @(#) ioapi library version 3.1 $".ljust(80);
    DI.EXEC_ID = "BCON_V5g_Darwin13_x86_64gfortran ".ljust(80);
    DI.FTYPE = 1 ;
    DI.CDATE = 2014031 ;
    DI.CTIME = 230334 ;
    DI.WDATE = 2014031 ;
    DI.WTIME = 230334 ;
    DI.SDATE = np.int32(metcroprops['SDATE']) ;
    DI.STIME = np.int32(metcroprops['STIME']) ;
    DI.TSTEP = np.int32(metcroprops['TSTEP']) ;
    DI.NTHIK = 1 ;
    DI.NCOLS = np.int32(metcroprops['NCOLS']) ;
    DI.NROWS = np.int32(metcroprops['NROWS']) ;
    DI.NLAYS = np.int32(metcroprops['NLAYS']) ;
    DI.NVARS = 106 ;
    DI.GDTYP = np.int32(metcroprops['GDTYP']) ;
    DI.P_ALP = np.float64(metcroprops['P_ALP']) ;
    DI.P_BET = np.float64(metcroprops['P_BET']) ;
    DI.P_GAM = np.float64(metcroprops['P_GAM']) ;
    DI.XCENT = np.float64(metcroprops['XCENT']) ;
    DI.YCENT = np.float64(metcroprops['YCENT']) ;
    DI.XORIG = np.float64(metcroprops['XORIG']) ;
    DI.YORIG = np.float64(metcroprops['YORIG']) ;
    DI.XCELL = np.float64(metcroprops['XCELL']) ;
    DI.YCELL = np.float64(metcroprops['YCELL']) ;
    DI.VGTYP = 7 ;
    DI.VGTOP = np.float32(metcroprops['VGTOP']) ;
    DI.VGLVLS = np.asarray(metcroprops['VGLVLS'], dtype = np.float32) ;
    DI.GDNAM = metcroprops['GDNAM'].ljust(16) ;
    DI.UPNAM = metcroprops['UPNAM'].ljust(16) ;
    setattr(DI, 'VAR-LIST', ''.join([i.ljust(16) for i in varlist + ["AIRDEN"]]));
    DI.FILEDESC = "ICON output file IC_CONC_1".ljust(80)
    DI.HISTORY = "" ;
    DI.sync()
    return DI

def makedefaultbcon(metbdyprops, args):
    """
    Create an empty bcon file whose dimensions
    are consistent with metbdyprops
    
    metbdyprops - dictionary of properties from IOAPI metadata
    """
    profilepath = args.BCONPROFILEPATH
    SPCLISTS = getspclists(args)
    GCLIST = SPCLISTS['GCLIST']
    AELIST = SPCLISTS['AELIST']
    SFLIST = SPCLISTS['SFLIST']
    NMLIST = SPCLISTS['NMLIST']
    if not profilepath is None:
        from PseudoNetCDF.cmaqfiles import bcon_profile
        bconfile = bcon_profile(profilepath)
        if not args.keepall: bconfile = getvarpnc(bconfile, GCLIST + AELIST + SFLIST + NMLIST)
        varlist = [varkey for varkey in bconfile.variables.keys() if varkey in (GCLIST + AELIST + SFLIST + NMLIST)]
        missingunits = [varkey for varkey, var in bconfile.variables.items() if var[:].ndim == 2 and varkey not in GCLIST + AELIST + NMLIST + SFLIST]
        missing = ', '.join(missingunits)
        if missing != '': warn(profilepath + ' contains ' + missing + ' whose units are assumed ppmV')
        GCLIST = missingunits + GCLIST
        def getvals(varkey):
            if varkey not in bconfile.variables:
                warn(varkey + ' not found in bcon profile path')
                return 1e-32
            invglvls = np.convolve([.5, .5], metbdyprops['VGLVLS'], mode = 'valid')
            outvglvls = np.convolve([.5, .5], bconfile.VGLVLS, mode = 'valid')
            invar = bconfile.variables[varkey]
            south = interpbound(invglvls, outvglvls, invar[:, 0])[:, None].repeat(metbdyprops['NCOLS'] + 1, 1)
            north = interpbound(invglvls, outvglvls, invar[:, 2])[:, None].repeat(metbdyprops['NCOLS'] + 1, 1)
            east = interpbound(invglvls, outvglvls, invar[:, 1])[:, None].repeat(metbdyprops['NROWS'] + 1, 1)
            west = interpbound(invglvls, outvglvls, invar[:, 3])[:, None].repeat(metbdyprops['NROWS'] + 1, 1)
            return np.concatenate([south, east, north, west], axis = 1)
    else:
        def getvals(varkey):
            return 1e-32
    DB = Dataset('dummybcon.nc', 'w', format = 'NETCDF3_CLASSIC')
    DB.createDimension('TSTEP', None)
    DB.createDimension('DATE-TIME', 2) ;
    DB.createDimension('VAR', 78) ;
    DB.createDimension('LAY', metbdyprops['NLAYS']) ;
    DB.createDimension('PERIM', metbdyprops['PERIM']) ;

    TFLAG = DB.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME')) ;
    TFLAG.units = "<YYYYDDD,HHMMSS>" ;
    TFLAG.long_name = "TFLAG           " ;
    TFLAG.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
    def getunit(varkey):
        if varkey in GCLIST:
            return 'ppmV'
        elif varkey in AELIST:
            return 'micrograms/m**3'
        elif varkey in NMLIST:
            return '#/m**3'
        elif varkey in SFLIST:
            return 'm**2/m**3'
    #varkey_unit = [(i, 'ppmV') for i in GCLIST] + [(i, 'micrograms/m**3') for i in AELIST] + [(i, '#/m**3') for i in NMLIST] + [(i, 'm**2/m**3') for i in SFLIST]
    varkey_unit = [(varkey, getunit(varkey)) for varkey in varlist]
    for varkey, unit in varkey_unit:
        var = DB.createVariable(varkey, 'f', ('TSTEP', 'LAY', 'PERIM'))
        var.long_name = varkey.ljust(16);
        var.units = unit.ljust(16);
        var.coordinates = "lon lat" ;
        var.var_desc = ("Variable %s" % varkey).ljust(80);
    
    for varkey, unit in varkey_unit:
        DB.variables[varkey][0] = getvals(varkey);
    
    TFLAG[0] = 0
    DB.IOAPI_VERSION = "$Id: @(#) ioapi library version 3.1 $".ljust(80);
    DB.EXEC_ID = "BCON_V5g_Darwin13_x86_64gfortran ".ljust(80);
    DB.FTYPE = 2 ;
    DB.CDATE = 2014031 ;
    DB.CTIME = 230334 ;
    DB.WDATE = 2014031 ;
    DB.WTIME = 230334 ;
    DB.SDATE = np.int32(metbdyprops['SDATE']) ;
    DB.STIME = np.int32(metbdyprops['STIME']) ;
    DB.TSTEP = np.int32(metbdyprops['TSTEP']) ;
    DB.NTHIK = 1 ;
    DB.NCOLS = np.int32(metbdyprops['NCOLS']) ;
    DB.NROWS = np.int32(metbdyprops['NROWS']) ;
    DB.NLAYS = np.int32(metbdyprops['NLAYS']) ;
    DB.NVARS = 78 ;
    DB.GDTYP = np.int32(metbdyprops['GDTYP']) ;
    DB.P_ALP = np.float64(metbdyprops['P_ALP']) ;
    DB.P_BET = np.float64(metbdyprops['P_BET']) ;
    DB.P_GAM = np.float64(metbdyprops['P_GAM']) ;
    DB.XCENT = np.float64(metbdyprops['XCENT']) ;
    DB.YCENT = np.float64(metbdyprops['YCENT']) ;
    DB.XORIG = np.float64(metbdyprops['XORIG']) ;
    DB.YORIG = np.float64(metbdyprops['YORIG']) ;
    DB.XCELL = np.float64(metbdyprops['XCELL']) ;
    DB.YCELL = np.float64(metbdyprops['YCELL']) ;
    DB.VGTYP = 7 ;
    DB.VGTOP = np.float32(metbdyprops['VGTOP']) ;
    DB.VGLVLS = np.asarray(metbdyprops['VGLVLS'], dtype = np.float32) ;
    DB.GDNAM = metbdyprops['GDNAM'].ljust(16) ;
    DB.UPNAM = metbdyprops['UPNAM'].ljust(16) ;
    setattr(DB, 'VAR-LIST', ''.join([i.ljust(16) for i in varlist + ["AIRDEN"]]));
    DB.FILEDESC = "BCON output file BNDY_CONC_1".ljust(80)
    DB.HISTORY = "" ;
    DB.sync()
    return DB


def get_template(option, gcversion):
    """
    Returns json template for pncglobal2cmaq
    
    Requries:
        option - string with mechanism name (only supports cb05) and CMAQ aerosol option (AE5 or AE6)
        gcversion - GEOS-Chem version (8, 9, 10, etc)
    """
    # First add the preamble
    out =r"""{
    "comment1": "comments are added like this",
    "comment2": "Each line has the form \"Species\": {\"expression\": \"some expression\", \"unit\": \"some unit\"},",
    "comment3": "Variables with the output unit (micrograms/m**3) will be multiplied by air density (moles/m**3) and molar mass in the script. It should not be done in the expression.",
    "comment4": "Set \"manual_unit\": true to bypass script unit corrections",
    "comment5": "If you do not have PSURF and/or TMPU, you can use values from the standard atmosphere by adding --expr=\"PSURF=np.ones_like(O3[:,0][:, None])*1013.25;TMPU=np.ones_like(O3[:])*np.array([287.7,286.9,286.0,285.2,284.3,283.4,282.6,281.7,280.7,279.8,278.9,277.9,276.8,275.3,273.6,271.8,270.0,268.2,265.8,262.8,259.7,256.4,252.9,249.2,245.2,241.0,236.4,230.5,223.6,216.8,216.6,216.6,216.6,216.6,216.6,216.6,216.6,217.5,219.8,222.0,225.3,233.6,250.5,269.2,260.0,237.7,214.3])[None, :47, None, None]\"",
    "comment6": "more comments are added like this",
    "PRESS": {
            "expression": "(hyam[:].reshape(1, -1).T + hybm[:].reshape(1, -1).T * PSURF[:][:, [0]].T).T * 100.",
            "outunit": "Pa"
        },
    "AIRMOLDEN": {
            "expression": "(hyam[:].reshape(1, -1).T + hybm[:].reshape(1, -1).T * PSURF[:][:, [0]].T).T * 100. / 8.3144621 / TMPU[:]",
            "outunit": "moles/m**3",
            
        }, 
    "AIRMASSDEN": {
            "expression": "0.0289645 * (hyam[:].reshape(1, -1).T + hybm[:].reshape(1, -1).T * PSURF[:][:, [0]].T).T * 100 / 8.3144621 / TMPU[:]",
            "outunit": "kg/m**3",
        },
    "CMAQSPECIES": {
        "AIRDEN": {
            "expression": "(hyam[:].reshape(1, -1).T + hybm[:].reshape(1, -1).T * PSURF[:][:, [0]].T).T * 100 / 8.3144621 / TMPU[:] * 6.022e23 / 1e6",
            "outunit": "molec/cm3",
            "manual_unit": true,
            "comment": "Full calculated unit",
        },"""
    # Next add the gas-phase
    issaprc07 = 'saprc07' in option
    iscb05 = 'cb05' in option
    iscb6  = 'cb6'  in option
    if iscb05 or iscb6:
        out += _CBSPCS
        if iscb05:
            out += _CB05SPCS 
        elif iscb6:
            warn('Check CB6 version and update nitrates accordingly')
            out += _CB6SPCS
            
        if 'ae6' in option:
            out += _GC_TO_AE6_COMMON
            if gcversion > 8:
                out += _GCgt8_to_AE6
            else:
                out += _GC8_to_AE6
    if out[-1] == ',':
        out = out[:-1]
    out += r"""
    }
}
"""
    
    return out

class speciesstruct(object):
    def __init__(self, name, ind, mode, density, version, found):
        """
        name - text string identifying aerosol
        ind - variable number (still needed?)
        mode - 0: Aitken, 1: Accumulation, 2: Coarse
        density - aerosol density (kg/m^3)
        version - 0: both, 5: AERO5, 6: AERO6
        """
        self.name = name
        self.mode = mode
        self.density = density
        self.version = version
        self.found = found

def repair_ae(f, myioo):
    """
    Create metavariables that are consistent with 
    mass from true variables.
    """
    if myioo is None:
        status = warn = error = warnings.warn
    else:
        warn = myioo.warn
        status = myioo.status
        error = myioo.error
    
    RHOSO4  = 1.8e3 # bulk density of aerosol sulfate
    RHONH4  = 1.8e3 # bulk density of aerosol ammonium
    RHONO3  = 1.8e3 # bulk density of aerosol nitrate
    RHOORG  = 1.3e3 # bulk density for aerosol organics following Carlton et al. 2010
    RHOSOIL = 2.6e3 # bulk density for aerosol soil dust
    RHOSEAS = 2.2e3 # bulk density for marine aerosol
    RHOANTH = 2.2e3 # bulk density for anthropogenic aerosol
    SGINIAT = 1.7   # initial sigma-G for Aitken mode
    SGINIAC = 2.0   # initial sigma-G for accumulation mode
    SGINICO = 2.2   # initial sigma-G for coarse mode
    DGINIAT = 0.01E-6  # geometric mean diameter for Aitken mode [ m ]
    DGINIAC = 0.07E-6  # geometric mean diameter for accum  mode [ m ]
    DGINICO = 1.0E-6   # geometric mean diameter for coarse mode [ m ]
    CONMIN  = 1.0E-30  # minimum concentration [ ug/m**3 ]
    nspecies = 57      # number of aerosol species treated

    #...conversion factors for number and surface area
    NUMFAC = {}
    NUMFAC['ATKN'] = 1.0 / ( ( DGINIAT ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINIAT ) ) ** 2.0 ) ) )
    NUMFAC['ACC']  = 1.0 / ( ( DGINIAC ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINIAC ) ) ** 2.0 ) ) )
    NUMFAC['COR']  = 1.0 / ( ( DGINICO ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINICO ) ) ** 2.0 ) ) )
    SRFFAC = {}
    SRFFAC['ATKN'] = pi / ( DGINIAT * exp( ( 5.0 / 2.0 ) * ( ( log( SGINIAT ) ) ** 2.0 ) ) )
    SRFFAC['ACC']  = pi / ( DGINIAC * exp( ( 5.0 / 2.0 ) * ( ( log( SGINIAC ) ) ** 2.0 ) ) )
    SRFFAC['COR']  = pi / ( DGINICO * exp( ( 5.0 / 2.0 ) * ( ( log( SGINICO ) ) ** 2.0 ) ) )

    bcspcs = [speciesstruct (    'ACLI',   0, 'ATKN', RHOSEAS, 0, False),
             speciesstruct (    'AECI',   0, 'ATKN', RHOANTH, 0, False),
             speciesstruct (    'ANAI',   0, 'ATKN', RHOSEAS, 0, False),
             speciesstruct (   'ANH4I',   0, 'ATKN',  RHONH4, 0, False),
             speciesstruct (   'ANO3I',   0, 'ATKN',  RHONO3, 0, False),
             speciesstruct (   'ASO4I',   0, 'ATKN',  RHOSO4, 0, False),
             speciesstruct (    'A25I',   0, 'ATKN', RHOANTH, 5, False),
             speciesstruct ( 'AORGPAI',   0, 'ATKN',  RHOORG, 5, False),
             speciesstruct (  'AOTHRI',   0, 'ATKN', RHOANTH, 6, False),
             speciesstruct ( 'APNCOMI',   0, 'ATKN',  RHOORG, 6, False),
             speciesstruct (   'APOCI',   0, 'ATKN',  RHOORG, 6, False),
             speciesstruct (   'AALKJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ABNZ1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ABNZ2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ABNZ3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (    'ACLJ',   0,  'ACC', RHOSEAS, 0, False),
             speciesstruct (    'AECJ',   0,  'ACC', RHOANTH, 0, False),
             speciesstruct (  'AISO1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AISO2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AISO3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (    'ANAJ',   0,  'ACC', RHOSEAS, 0, False),
             speciesstruct (   'ANH4J',   0,  'ACC',  RHONH4, 0, False),
             speciesstruct (   'ANO3J',   0,  'ACC',  RHONO3, 0, False),
             speciesstruct (  'AOLGAJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AOLGBJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AORGCJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (   'ASO4J',   0,  'ACC',  RHOSO4, 0, False),
             speciesstruct (   'ASQTJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATOL1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATOL2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATOL3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATRP1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATRP2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AXYL1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AXYL2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AXYL3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (    'A25J',   0,  'ACC', RHOANTH, 5, False),
             speciesstruct ( 'AORGPAJ',   0,  'ACC',  RHOORG, 5, False),
             speciesstruct (    'AALJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'ACAJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'AFEJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (     'AKJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'AMGJ',   0,  'ACC', RHOSEAS, 6, False),
             speciesstruct (    'AMNJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (  'AOTHRJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct ( 'APNCOMJ',   0,  'ACC',  RHOORG, 6, False),
             speciesstruct (   'APOCJ',   0,  'ACC',  RHOORG, 6, False),
             speciesstruct (    'ASIJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'ATIJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'ACLK',   0,  'COR', RHOSEAS, 0, False),
             speciesstruct (   'ACORS',   0,  'COR', RHOANTH, 0, False),
             speciesstruct (   'ANH4K',   0,  'COR',  RHONH4, 0, False),
             speciesstruct (   'ANO3K',   0,  'COR',  RHONO3, 0, False),
             speciesstruct (   'ASO4K',   0,  'COR',  RHOSO4, 0, False),
             speciesstruct (   'ASOIL',   0,  'COR', RHOSOIL, 0, False),
             speciesstruct (    'ANAK',   0,  'COR', RHOSEAS, 5, False),
             speciesstruct ( 'ASEACAT',   0,  'COR', RHOSEAS, 6, False)]

    for spc in bcspcs:
        try:
            data_shape = f.variables[spc.name].shape
            break
        except:
            pass
    else:
        warn("There are no aerosol species")
        return
    not_found = defaultdict(lambda: [])
    for spc in bcspcs:
        try:
            spc.found = spc.name in f.variables.keys()
        except:
            spc.found = False
        if not spc.found:
            not_found[spc.version].append(spc.name)
    version_check = []
    for k, v in not_found.items():
        if k != 0:
            version_check.append(k)
        warn('Some variables from %d were not found: %s' % (k, v))
    
    if len(version_check) > 1:
        warn('Some variables from aerosol versions %s were not found' % ' and '.join([str(v) for v in version_check]))
    moment3 = dict([(k, zeros(data_shape, dtype = 'f')) for k in 'ATKN ACC COR'.split()])
    
    for spc in bcspcs:
        try:
            bcval = f.variables[spc.name]
        except KeyError:
            continue
        v = 1.0e-9*6.0/( pi*spc.density ) * bcval[:]
        moment3[spc.mode] += v
    
    for modek, modv in moment3.items():
        numkey = 'NUM' + modek
        srfkey = 'SRF' + modek
        if numkey in f.variables.keys():
            numvar = f.variables[numkey]
        else:
            bcspc = [spc.name for spc in bcspcs if spc.found][0]
            dims = f.variables[bcspc].dimensions
            numvar = f.createVariable(numkey, 'f', dims)
            numvar.units = '#/m**3'.ljust(16);
            numvar.long_name = numkey.ljust(16);
            numvar.var_desc = numkey.ljust(16);

        if srfkey in f.variables.keys():
            srfvar = f.variables[srfkey]
        else:
            bcspc = [spc.name for spc in bcspcs if spc.found][0]
            dims = f.variables[bcspc].dimensions
            srfvar = f.createVariable(srfkey, 'f', dims)
            srfvar.units = 'm**2/m**3'.ljust(16);
            srfvar.long_name = srfkey.ljust(16);
            srfvar.var_desc = srfkey.ljust(16);
        numvar[:] = NUMFAC[modek] * modv
        srfvar[:] = SRFFAC[modek] * modv
    f.sync()

def makeregriddedbdy(metbdy, args, spcs):
    """
    Requires:
        metbdy - METBDY3D output from MCIP with longitude and latitude variables
        args - args from script argparse
        spcs - list of species
    Outputs:
        list netcdf-like files with ND49 data regridded to metbdy domain
    """
    lonlatcoords = '/'.join(['%f,%f' % (lon, lat) for lon, lat in zip(metbdy.variables['longitude'][:].ravel(), metbdy.variables['latitude'][:].ravel())])
    out = []
    for ND49, ND49_REGRID_BDY in zip(args.ND49, args.ND49_REGRID_BDY):
        outf = extract(getvarpnc(ND49, None), lonlatcoords, method = args.extractmethod)
        if args.persistintermediate:
            pncgen(outf, ND49_REGRID_BDY, verbose = False)
            outf = Dataset(ND49_REGRID_BDY)
        out.append(outf)
    return out

def makeregriddedcro(metcro, args, spcs):
    """
    Requires:
        metcro - METCRO3D output from MCIP with longitude and latitude variables
        args - args from script argparse
        spcs - list of species
    Outputs:
        otu - list netcdf-like files with ND49 data regridded to metcro domain
    """
    lonlatcoords = '/'.join(['%f,%f' % (lon, lat) for lon, lat in zip(metcro.variables['longitude'][:].ravel(), metcro.variables['latitude'][:].ravel())])
    out = []
    #import pdb; pdb.set_trace()
    for ND49, ND49_REGRID_CRO in zip(args.ND49, args.ND49_REGRID_CRO):
        outf = extract(slice_dim(getvarpnc(ND49, None), 'time,0'), lonlatcoords, method = args.extractmethod)
        if args.persistintermediate:
            pncgen(outf, ND49_REGRID_CRO, verbose = False)
            outf = Dataset(ND49_REGRID_CRO)
        out.append(outf)
        
    return out

def getdefault(oldcon, vark, noutstep):
    """
    Requires:
        oldcon - netcdf-like file with default values for [IB]CON file
        vark - variable name (must be in oldcon.variables)
        noutstep - number of output steps
    Output:
        defval - array of default values for vark
    """
    defval = np.ma.filled(oldcon.variables[vark][:], 0)
    if defval.shape[0] == 1:
        defval = defval.repeat(noutstep, 0)
    elif defval.shape[0] > noutstep:
        warn('Default (I,B)CON has time dimension that is greater than output (%d>%d); only first %d are used.' % (defval.shape[0], noutstep, noutstep))
        defval = defval[0:noutstep]
    return defval


def makeibcon(args):
    """
    Main logic to make [IB]CON files
    """
    global messages
    if args.verbose:
        print('Starting Main')
    messages += ' '.join(sys.argv[:]) + '\n'
    mappings_file = json.load(open(args.mapping, mode = 'r'))
    mappings = mappings_file['CMAQSPECIES']
    spcs_list = []
    for spcs in [compile(mapping['expression'], 'mapping', 'eval').co_names for mapping in [mv for mv in mappings.values()] + [mappings_file['AIRMOLDEN']]]:
        spcs_list.extend(spcs)
    
    spcs = ','.join(spcs_list)
    if args.METBDY3D is None:
        dobcon = False
    else:
        dobcon = True

    if args.METCRO3D is None:
        doicon = False
    else:
        doicon = True
    
    if dobcon:
        if args.verbose:
            print('Starting BCON Prep')
        metbdyfiles, metbdyargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.METBDY3D.split(' '), parser = None)
        try: metbdy = getvarpnc(metbdyfiles[0], ['PRES'])
        except: metbdy = getvarpnc(metbdyfiles[0], None)
        add_cf_from_ioapi(metbdy)
        regridded_nd49_bdy = makeregriddedbdy(metbdy, args, spcs)
        

        metbdyprops = dict([(propk, getattr(metbdy, propk)) for propk in metbdy.ncattrs()])
        metbdyprops['VGLVLSCDL'] = 'f, '.join(['%f' % vgl for vgl in metbdyprops['VGLVLS']]) + 'f'
        metbdyprops['PERIM'] = len(metbdy.dimensions['PERIM'])
        metbdyprops['GDTYP'] = metbdy.GDTYP

        if args.BCON == 'dummybcon.nc':
            oldbcon = makedefaultbcon(metbdyprops, args)
        else:
            bconfiles, bconargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.BCON.split(' '), parser = None)
            oldbcon = bconfiles[0]
        newbcon = Dataset(args.NEWBCON, mode = 'ws', format = 'NETCDF3_CLASSIC')
        for propk in oldbcon.ncattrs():
            propv = getattr(oldbcon, propk)
            setattr(newbcon, propk, propv)

        newbcon.createDimension('TSTEP', None)
        newbcon.createDimension('DATE-TIME', 2)
        for dimk in ('LAY', 'PERIM'):
            newbcon.createDimension(dimk, len(oldbcon.dimensions[dimk]))
        for vark, oldv in oldbcon.variables.items():
            if vark == 'TFLAG': continue
            newv = newbcon.createVariable(vark, oldv.dtype.char, oldv.dimensions)
            for propk in oldv.ncattrs():
                propv = getattr(oldv, propk)
                setattr(newv, propk, propv)
        try:
            nbconoutsteps = sum([len(tmpf.dimensions['time']) for tmpf in regridded_nd49_bdy])
        except:
            nbconoutsteps = sum([len(tmpf.dimensions['TSTEP']) for tmpf in regridded_nd49_bdy])
        if len(newbcon.dimensions['PERIM']) != metbdyprops['PERIM']:
            raise ValueError('Default (I,B)CON has different perimeter dimension (%d) than the output file (%d).' % (len(newbcon.dimensions['PERIM']), metbdyprops['PERIM']))
        obconsteps = len(oldbcon.dimensions['TSTEP'])
        if obconsteps != 1 and obconsteps < nbconoutsteps:
            raise ValueError('Default BCON has time dimension that is less than output (%d<%d) and not 1. Use time-independent file or file with same times.' % (obconsteps, nbconoutsteps))

    
    if doicon:
        if args.verbose:
            print('Starting ICON Prep')
        metcrofiles, metcroargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.METCRO3D.split(' '), parser = None)
        try: metcro = slice_dim(getvarpnc(metcrofiles[0], ['PRES']), 'TSTEP,0')
        except: metcro = slice_dim(getvarpnc(metcrofiles[0], None), 'TSTEP,0')
        add_cf_from_ioapi(metcro)
        regridded_nd49_cro = makeregriddedcro(metcro, args, spcs)
        metcroprops = dict([(propk, getattr(metcro, propk)) for propk in metcro.ncattrs()])
        metcroprops['VGLVLSCDL'] = 'f, '.join(['%f' % vgl for vgl in metcroprops['VGLVLS']]) + 'f'
        metcroprops['GDTYP'] = metcro.GDTYP

        if args.ICON == 'dummyicon.nc':
            oldicon = makedefaulticon(metcroprops, args)
        else:
            iconfiles, iconargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.ICON.split(' '), parser = None)
            oldicon = iconfiles[0]
        oldicon = Dataset(args.ICON, mode = 'r+s', format = 'NETCDF3_CLASSIC')
        newicon = Dataset(args.NEWICON, mode = 'ws', format = 'NETCDF3_CLASSIC')
        for propk in oldicon.ncattrs():
            propv = getattr(oldicon, propk)
            setattr(newicon, propk, propv)
    
        newicon.createDimension('TSTEP', None)
        for dimk in ('DATE-TIME', 'LAY', 'ROW', 'COL'):
            newicon.createDimension(dimk, len(oldicon.dimensions[dimk]))
        for vark, oldv in oldicon.variables.items():
            if vark == 'TFLAG': continue
            newv = newicon.createVariable(vark, oldv.dtype.char, oldv.dimensions)
            for propk in oldv.ncattrs():
                propv = getattr(oldv, propk)
                setattr(newv, propk, propv)
        if len(newicon.dimensions['ROW']) != metcroprops['NROWS']:
            raise ValueError('Default (I,B)CON has different ROW dimension (%d) than the output file (%d).' % (len(newicon.dimensions['ROW']), metcroprops['NROWS']))
        if len(newicon.dimensions['COL']) != metcroprops['NCOLS']:
            raise ValueError('Default (I,B)CON has different COL dimension (%d) than the output file (%d).' % (len(newicon.dimensions['COL']), metcroprops['NCOLS']))
        

    minval = 1e-32
    addedkeys = []
    infiles = []
    if dobcon:
        infiles += [(newbcon, ('TSTEP', 'LAY', 'PERIM'))]
    if doicon:
        infiles += [(newicon, ('TSTEP', 'LAY', 'ROW', 'COL'))]
    if args.verbose:
        print('Starting Var Def')
    for vark, varo in mappings.items():
        for newcon, dims in infiles:
            if vark not in newcon.variables.keys():
                newv = newcon.createVariable(vark, 'f', dims)
                newv.units = varo['outunit'].ljust(16)
                newv.var_desc = newv.long_name = vark.ljust(16)
                addedkeys.append(vark)

    infiles = []
    if dobcon:
        infiles += [(metbdy, regridded_nd49_bdy, oldbcon, newbcon, nbconoutsteps)]
    if doicon:
        infiles +=  [(metcro, regridded_nd49_cro, oldicon, newicon, 1)]
    general_messages = messages
    messages = ""
    for metfile, regridded_nd49, oldcon, newcon, noutstep in infiles:
        if args.verbose:
            print('Starting Var Calc')
        if args.sigmaeta or args.sigma:
            cpressv = np.convolve([0.5, 0.5], metfile.VGLVLS, mode = 'valid')
        else:
            cpressv = metfile.variables['PRES']
            assert(cpressv.units.strip() == 'Pa')
        oldkeys = []
        for vark, varo in newcon.variables.items():
            if args.verbose:
                print('Starting', vark, 'Calc')
            toff = 0
            if vark in ('TFLAG',): continue
            if vark in mappings:
                sys.stdout.write(vark + ', ')
                sys.stdout.flush()
                expr, ounit = eval('expression, outunit', None, mappings[vark])
                # Temporary version for units
                nd49 = regridded_nd49[0]
                coexpr = compile(expr, 'expr', 'eval')
                inunits = [nd49.variables[cn].units.strip() for cn in coexpr.co_names if cn in nd49.variables]
                if len(inunits) == 0:
                    warn("Cannot evaluate any parts of '%s' = '%s'; using default file" % (vark, expr))
                    oldkeys.append(vark)
                    minout = getdefault(oldcon, vark, noutstep)
                else:
                    inC = dict([(cn, getattr(nd49.variables[cn], 'carbon', 1.)) for cn in coexpr.co_names if cn in nd49.variables])
                    inkgpermole = dict([(cn, getattr(nd49.variables[cn], 'kgpermole', np.nan)) for cn in coexpr.co_names if cn in nd49.variables])
                    def conv(matcho):
                        found = {}
                        found.update(matcho.groupdict())
                        spc = found['spc']
                        if spc not in nd49.variables:
                            return spc
                        found['carbon'] = inC[spc]
                        found['kgpermole'] = inkgpermole[spc]
                        if 'micrograms' in ounit:
                            out =  '%(spc)s[:] / %(carbon)f * %(kgpermole)f' % found
                        elif found['carbon'] == 1.:
                            out =  '%(spc)s[:]' % found
                        else:
                            out =  '%(spc)s[:] / %(carbon)f' % found
                        return out
                    for cn in coexpr.co_names:
                        expr = re.compile(r'(?P<spc>\b' + cn + r'\b)').sub(conv, expr)
            
                    inunit = inunits[0]
                    manual_unit = mappings[vark].get('manual_unit', False)
                    ounit = ounit.strip()
                    if args.verbose:
                        sys.stdout.write('\n%s %s\n' % (expr, ounit))
                        sys.stdout.flush()
                    out = np.zeros((noutstep,) + varo[:].shape[1:], dtype = 'f')
                    for ndi, nd49 in enumerate(regridded_nd49):
                        if not args.sigmaeta:
                            rpressv = eval(mappings_file['PRESS']['expression'], None, nd49.variables)
                        
                        if not (ounit == varo.units.strip()):
                            sys.stdout.write('\n%s %s %s\n' % (vark, varo.units.strip(), ounit))
                            sys.stdout.flush()
                        
                        try:
                            temp_val = eval(expr, None, nd49.variables)[:]
                        except:
                            warn("Cannot evaluate part of '%s' = '%s'" % (vark, expr))
                            oldkeys.append(vark)
                            out = getdefault(oldcon, vark, noutstep)
                            break
                    
                        if manual_unit:
                            pass
                        else:
                            assert((np.array(inunits) == inunit).all())
                            if inunit == ounit:
                                pass
                            elif ounit == 'ppmV' and inunit in ('pptv', 'pptC'):
                                # pptC has automatically been converted to pptv
                                temp_val *= 1e-6
                            elif ounit == 'ppmV' and inunit in ('ppbv', 'ppbC'):
                                # ppbC has automatically been converted to ppbv
                                temp_val *= 1e-3
                            elif ounit == 'micrograms/m**3' and (inunit == 'ppbv' or inunit == 'pptv'):
                                airmolden = eval(mappings_file['AIRMOLDEN']['expression'], None, nd49.variables)
                                temp_val *= airmolden
                                if inunit == 'pptv':
                                    temp_val *= 1e-3
                            else:
                                raise ValueError('Error: in unit/outunit combo unknown: "%s", "%s"' % (inunit, ounit))
                        toff
                        oldrm = 0
                        for ti, temp_hour in enumerate(temp_val):
                            if args.verbose:
                                print('\nHour', ti, ':      ', end = '')
                            if args.sigmaeta or args.sigma:
                                cpress = cpressv[:, None].repeat(out[0,0].size, 1)
                            else:
                                cpress = cpressv[ti, :]
                                cpress = cpress.reshape(cpress.shape[0], -1)
                                assert((np.diff(cpress, axis = 0).mean(1) < 0).all())
                            if args.sigmaeta:
                                if hasattr(nd49, 'VGLVLS'):
                                    rpress = ((nd49.VGLVLS[:-1] + nd49.VGLVLS[1:]) * .5)[:, None].repeat(out[0, 0].size, 1)
                                else:
                                    rpress = ((nd49.variables['etam_pressure'] - newcon.VGTOP / 100) / (1013.25 - newcon.VGTOP / 100))[:, None].repeat(out[0,0].size, 1)
                            else:
                                if mappings_file['PRESS']['outunit'] == 'Pa':
                                    rpress = rpressv[ti,:]
                                elif mappings_file['PRESS']['outunit'] == 'hPa':
                                    rpress = rpressv[ti,:] * 100
                                else:
                                    warn("Assuming pressure is Pa, but got %s" % mappings_file['PSURF']['outunit'])
                                if args.sigma:
                                    rpress = (rpress - newcon.VGTOP) / (100 * nd49.variables['PSURF'][ti,0][None].reshape(1, -1) - newcon.VGTOP)
                            rpress = rpress.reshape(rpress.shape[0], cpress.shape[1])
                            assert((np.diff(rpress, axis = 0).mean(1) < 0).all())
                                                
                            outvals = np.zeros(cpress.shape, dtype = 'f')
                            for pi, (cp, rp, column) in enumerate(zip(cpress.T, rpress.T, temp_hour.T)):
                                outvals[:, pi] = np.interp(cp[::-1], rp[::-1], column[::-1])[::-1]
                                if args.verbose:
                                    print(4*'\b',end = '')
                                    print('%3.0f%%' % (pi/cpress.shape[1]*100),end = '')
                            out[toff + ti, :, :] = outvals.reshape(*out.shape[1:])
                        toff = toff + ti + 1
                    minout = np.maximum(out, minval)
                    if args.verbose:
                        sys.stdout.write('\nSPC     0    25    50    75   100')
                        sys.stdout.write('\n%s %.5e %.5e %.5e %.5e %.5e' % tuple(['Raw'] + np.percentile(out, [0, 25, 50, 75, 100]).tolist()))
                        sys.stdout.flush()
                        sys.stdout.write('\n%s %.5e %.5e %.5e %.5e %.5e\n' % tuple(['Min'] + np.percentile(minout, [0, 25, 50, 75, 100]).tolist()))
                        sys.stdout.flush()
            else:
                oldkeys.append(vark)
                minout = getdefault(oldcon, vark, noutstep)
            if args.verbose:
                print('Writing out', vark, ti)
            varo[0:minout.shape[0], :, :] = minout[:]

        #varlist = getattr(newcon, 'VAR-LIST')
        #nvars = getattr(newcon, 'NVARS') + len(addedkeys)
        #setattr(newcon, 'VAR-LIST', varlist + ''.join([nk.ljust(16) for nk in addedkeys]))
        allkeys = [nk.ljust(16) for nk in newcon.variables.keys() if nk != 'TFLAG']
        nvars = len(allkeys)
        newvarlist = ''.join(allkeys)
        setattr(newcon, 'VAR-LIST', newvarlist)
        setattr(newcon, 'NVARS', nvars)
        newcon.createDimension('VAR', nvars)
        tflag = newcon.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
        tflag.units = "<YYYYDDD,HHMMSS>" ;
        tflag.long_name = "TFLAG           " ;
        tflag.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
        if args.timeindependent:
            newcon.SDATE = np.int32(-635)
            newcon.STIME = np.int32(0)
            newcon.TSTEP = np.int32(0)
            tflag[:] = 0.
        else:
            gtime = gettimes(regridded_nd49[0])
            gsdate = gtime[0].strftime('%Y-%m-%d %H:%M:%S')
            if args.sdate is None:
                args.sdate = gsdate
            if args.tstep is None:
                try:
                    total_seconds = (gtime[1] - gtime[0]).total_seconds()
                    hours, remainder = divmod(total_seconds,60*60)
                    minutes, seconds = divmod(remainder,60)
                    args.tstep = '%02d%02d%02d' % (hours, minutes, seconds)
                except IndexError:
                    raise ValueError('Unable to determine timestep from file; use --tstep option to assign manually')
            sdate = datetime.strptime(args.sdate, '%Y-%m-%d %H:%M:%S')
            newcon.SDATE = np.int32(sdate.strftime('%Y%j'))
            newcon.STIME = np.int32(sdate.strftime('%H%M%S'))
            newcon.TSTEP = np.int32(args.tstep)
            sdate = datetime.strptime(str(newcon.SDATE) + ' ' + '%06d' % newcon.STIME, '%Y%j %H%M%S')
            dt = timedelta(hours = newcon.TSTEP // 10000.) + (datetime.strptime('%06d' % (newcon.TSTEP % 10000.), '%H%M%S') - datetime(1900, 1, 1))
            tnow = sdate
            for ti in range(varo.shape[0]):
                jdate = int(tnow.strftime('%Y%j'))
                itime = int(tnow.strftime('%H%M%S'))
                tflag[ti, :, :] = np.array([[jdate, itime]]).repeat(newcon.NVARS, 0)
                tnow = tnow + dt
        newcon.sync()
        repair_ae(newcon, None)
        if len(addedkeys) > 0:
            warn('CON did not originally have:\n' + ', '.join(addedkeys) + '\n')
        if len(oldkeys) > 0:
            warn('Using old CON for:\n' + ', '.join(oldkeys) + '\n')
        setattr(newcon, 'HISTORY', getattr(newcon, 'HISTORY', '') + general_messages + messages)
        newcon.sync()

from PseudoNetCDF.sci_var import stack_files, getvarpnc, PseudoNetCDFFile

class bpchwithstdatm(PseudoNetCDFFile):
    """
    Designed to add PSURF and TMPU variables based on simple standard atmospher assumption.
    """
    def __init__(self, bpath, vertgrid = 'GEOS-5-REDUCED', nogroup = True, verbose = False):
        PseudoNetCDFFile.__init__(self)
        from PseudoNetCDF.geoschemfiles import bpch
        self.variables = {}
        self.dimensions = {}
        infile = self._infile = bpch(bpath, vertgrid = vertgrid, nogroup = nogroup)
        NTIMES = len(infile.dimensions['time'])
        NLAYS = len(infile.dimensions['layer'])
        NLATS = len(infile.dimensions['latitude'])
        NLONS = len(infile.dimensions['longitude'])
        o = np.ones((NTIMES, NLAYS, NLATS, NLONS), dtype = 'f')
        if 'layer1' not in infile.dimensions:
            infile.createDimension('layer1', 1)
        pvar = infile.createVariable('PSURF', 'f', ('time', 'layer1', 'latitude', 'longitude'))
        pvar.units = 'hPa'
        pvar[:] = 1013.25
        if vertgrid == 'GEOS-5-REDUCED':
            layer = 'layer47'
        elif vertgrid == 'GEOS-5-NATIVE':
            layer = 'layer72'
        else:
            raise KeyError('vertgrid must be GEOS-5-NATIVE or GEOS-4-NATIVE')
        
        tvar = infile.createVariable('TMPU', 'f', ('time', 'layer', 'latitude', 'longitude'))
        tvar.units = 'K'
        tvar[:] = np.array([287.7,286.9,286.0,285.2,284.3,283.4,282.6,281.7,280.7,279.8,278.9,277.9,276.8,275.3,273.6,271.8,270.0,268.2,265.8,262.8,259.7,256.4,252.9,249.2,245.2,241.0,236.4,230.5,223.6,216.8,216.6,216.6,216.6,216.6,216.6,216.6,216.6,217.5,219.8,222.0,225.3,233.6,250.5,269.2,260.0,237.7,214.3])[None, :NLAYS, None, None]*o
        pncgen(infile, self, verbose = verbose)
    
        
        
class benchmark(PseudoNetCDFFile):
    """
    Designed to work with CTM runs with PEDGE-$_PSURF and DAO-3D-$_TMPU and IJ-AVG-$_*
    
    Available here:
    http://ftp.as.harvard.edu/gcgrid/geos-chem/1yr_benchmarks/
    """
    def __init__(self, bpath, vertgrid = 'GEOS-5-NATIVE'):
        PseudoNetCDFFile.__init__(self)
        from PseudoNetCDF.geoschemfiles import bpch
        inf = bpch(bpath, vertgrid = vertgrid)
        group_varkeys = [('', 'time'), ('', 'latitude'), ('', 'longitude'), ('', 'hyam'), ('', 'hybm'), ('', 'etam_pressure'), ('PEDGE-$', 'PSURF'), ('DAO-3D-$', 'TMPU')]
        group_varkeys += [('IJ-AVG-$', k.replace('IJ-AVG-$_', '')) for k in inf.variables.keys() if 'IJ-AVG-$' in k]
        getvarpnc_str = [((group + '_' + varkey) if group != '' else varkey) for group, varkey in group_varkeys]
        inf = getvarpnc(inf, getvarpnc_str)
        outf = self
        outf.createDimension('time', len(inf.dimensions['time']))
        nlay = len(inf.dimensions['layer'])
        outf.createDimension('layer', nlay)
        outf.createDimension('layer%d' % nlay, nlay)
        outf.createDimension('layer%d' % (nlay + 1), nlay + 1)
        outf.createDimension('latitude', len(inf.dimensions['latitude']))
        outf.createDimension('longitude', len(inf.dimensions['longitude']))
        outf.createDimension('nv', 2)
        outf.createDimension('tnv', 2)
        
        for group, varkey in group_varkeys:
            inkey = ('' if group == '' else (group + '_')) + varkey
            invar = inf.variables[inkey]
            outvar = outf.createVariable(varkey, invar.dtype.char, invar.dimensions)
            for pk in invar.ncattrs():
                setattr(outvar, pk, getattr(invar, pk))
            outvar[:] = invar[:]
        
        outf.sync()
registerreader('benchmark', benchmark)
registerreader('bpchwithstdatm', bpchwithstdatm)

if __name__ == '__main__':
    from PseudoNetCDF.pncparse import getparser, pncparse
    
    parser = getparser(has_ofile = False, plot_options = False, interactive = False)
    parser.add_argument('--withnco', action = 'store_true', dest = 'withnco', default = False, help = 'Use NCO to make dummy BCON/ICON.')
    parser.add_argument('--BCON', default = "dummybcon.nc", type = str, help='path (or PseudoNetCDF commands) to an old BCON file (default dummybcon.nc - created on the fly)')
    parser.add_argument('--ICON', default = "dummyicon.nc", help='path (or PseudoNetCDF commands) to an old ICON file')
    parser.add_argument('--persist', dest = 'persistintermediate', action = 'store_true', default = False, help = 'Save interpolated files (may have speed benefits)')
    parser.add_argument('--sigmaeta', action = 'store_true', help = 'Use sigma from CMAQ and calculate sigma for GEOS-Chem from pure eta levels for interpolation')
    parser.add_argument('--sigma', action = 'store_true', help = 'Use sigma from CMAQ and calculate sigma for GEOS-Chem from pressure levels for interpolation')
    parser.add_argument('--timeindependent', default = False, action = 'store_true', help = 'Start date/time is any start date time')
    parser.add_argument('-d', '--sdate', default = None, help = 'Start date in YYYY-MM-DD HH:MM:SS format')
    parser.add_argument('--METBDY3D', default = None, type = str, help='path (or PseudoNetCDF commands) to a MCIP METBDY3D file')
    parser.add_argument('--METCRO3D', default = None, help='path (or PseudoNetCDF commands) to a MCIP METCRO3D file')
    parser.add_argument('--AENML', default = None, help='Aerosol namelist')
    parser.add_argument('--NRNML', default = None, help='Nonreactive namelist')
    parser.add_argument('--GCNML', default = None, help='Gas concentration namelist')
    parser.add_argument('--keepall', default = False, action = 'store_true', help = 'Keep all species from profile data, even if not explicitly used')
    parser.add_argument('--BCONPROFILEPATH', default = None, help='path to BCON profile.dat')
    parser.add_argument('--ICONPROFILEPATH', default = None, help='path to ICON profile.dat')
    parser.add_argument('--tstep', default = None, help = 'Time increment between ND49 files')
    parser.add_argument('--mapping', default = 'mappings.json', help = 'Path to mappings file (i.e., json formatted dictionary); use --template to get a mapping')
    parser.add_argument('--outbfolder', default = ".", help = 'Path to output folder for BCON.')
    parser.add_argument('--outifolder', default = ".", help = 'Path to output folder for ICON.')
    parser.add_argument('--template', dest = 'template', default = None, choices = ['cb6mp_ae6_aq', 'cb05tucl_ae6_aq', 'cb05tucl_ae5_aq', 'saprc07tc_ae6_aq', 'saprc07tc_ae5_aq'], type = str, help = 'Print template mappings.json to stdout and exit (still requires ND49 argument)')
    parser.add_argument('--gc-version', dest = 'gcversion', default = 9, choices = [8, 9, 10, 11], type = int, help = 'GEOS-Chem version is used to determine SOA mapping. Above 8 (9-11) use the same SOA species.')
    parser.set_defaults(format = 'bpch,nogroup=True')
    parser.description = """pncgeos2cmaq makes boundary conditions from GEOS-Chem for CMAQ."""
    parser.epilog = """
Requirements:
    ifile - path to either a ND49 or TPCORE
    mapping - path to a table of mappings (see --template for defaults)
    one or both of:
        --METBDY3D - path to METBDY3D MCIP output file
        --METCRO3D - path to METCRO3D MCIP output file

Example:
    $ pncglobal2cmaq.py inputs/ts20120301.bpch --template=cb05tucl_ae6_aq > mappings.json
    $ pncglobal2cmaq.py inputs/ts20120301.bpch --METBDY3D inputs/METBDY3D_20120301 --METCRO3D inputs/METCRO3D_20120301 --mapping=mappings.json  --outifolder=outputs/ --outbfolder=outputs/
    
Example where BPCH file is missing PSURF and TMPU
    
"""
    
    args = parser.parse_args()
    if not args.template is None:
        print(get_template(args.template, args.gcversion))
        exit()
    if not os.path.exists(args.mapping):
        parser.print_help()
        raise ValueError("mappings.json must exist and can be created using the --template option")
    ifiles, args = pncparse(has_ofile = False, parser = parser)
    args.ND49 = ifiles
    args.ND49PATH = args.ipath
    if args.METCRO3D is None and args.METBDY3D is None:
        parser.print_help()
        raise ValueError('METCRO3D or METBDY3D are required (or both)')
    args.ND49_REGRID_BDY = [os.path.join(args.outbfolder, os.path.basename(p).replace('.nc', '') + '.BDY.nc') for p in args.ND49PATH]
    args.ND49_REGRID_CRO = [os.path.join(args.outifolder, os.path.basename(p).replace('.nc', '') + '.CRO.nc') for p in args.ND49PATH]
    args.NEWBCON = os.path.join(args.outbfolder, os.path.basename(args.ND49PATH[0].replace('.nc', '') + '.BCON.nc'))
    args.NEWICON = os.path.join(args.outifolder, os.path.basename(args.ND49PATH[0].replace('.nc', '') + '.ICON.nc'))
    makeibcon(args)
