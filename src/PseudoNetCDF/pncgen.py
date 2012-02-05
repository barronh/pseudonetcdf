from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(f,outpath, inmode = 'r', outmode = 'w'):
    p2n = Pseudo2NetCDF()
    return p2n.convert(f, outpath, inmode = inmode, outmode = outmode)
    
