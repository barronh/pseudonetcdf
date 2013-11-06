from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(f,outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4_CLASSIC'):
    p2n = Pseudo2NetCDF()
    return p2n.convert(f, outpath, inmode = inmode, outmode = outmode, format = format)
    
if __name__ == '__main__':
    from pncparse import pncparser
    ifile, ofile, options = pncparser()
    pncgen(f, ofile, outmode = options.mode, format = options.outformat)
