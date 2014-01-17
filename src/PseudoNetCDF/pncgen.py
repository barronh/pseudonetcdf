from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(ifile,outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4_CLASSIC', verbose = True):
    p2n = Pseudo2NetCDF()
    p2n.verbose = verbose
    return p2n.convert(ifile, outpath, inmode = inmode, outmode = outmode, format = format)
    
def main():
    from pncparse import pncparser
    ifile, ofile, options = pncparser(has_ofile = True)
    pncgen(ifile, ofile, outmode = options.mode, format = options.outformat, verbose = options.verbose)

if __name__ == '__main__':
    main()
