from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(ifile,outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4_CLASSIC', verbose = True):
    p2n = Pseudo2NetCDF()
    p2n.verbose = verbose
    return p2n.convert(ifile, outpath, inmode = inmode, outmode = outmode, format = format)
    
if __name__ == '__main__':
    from pncparse import pncparser
    ifile, ofile, options = pncparser()
    pncgen(ifile, ofile, outmode = options.mode, format = options.outformat, verbose = options.verbose)
