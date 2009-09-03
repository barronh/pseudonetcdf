from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(f,outpath, mode='w'):
    p2n = Pseudo2NetCDF()
    p2n.convert(f, outpath)
    