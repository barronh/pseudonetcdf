from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(f,outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4'):
    p2n = Pseudo2NetCDF()
    return p2n.convert(f, outpath, inmode = inmode, outmode = outmode, format = format)
    
if __name__ == '__main__':
    from optparse import OptionParser
    from camxfiles.Memmaps import uamiv
    from icarttfiles.ffi1001 import ffi1001
    from geoschemfiles import bpch

    parser = OptionParser()
    parser.set_usage("""Usage: python -m %prog [-f uamiv|bpch|ffi1001] ifile [ofile]

    ifile - path to a file formatted as type -f
    ofile - path to the desired output
    
    -f --format - format of the file either uamiv (CAMx), bpch (GEOS-Chem) or ffi1001 (optionally has comma delimited arguments for opening the format)
    """)

    parser.add_option("-f", "--format", dest = "format", default = 'uamiv', help = "File format")
    (options, args) = parser.parse_args()
    
    if len(args) == 0 or len(args) > 2:
        parser.print_help()
        exit()
    
    ifile = args[0]
    if len(args) == 1:
        ofile = ifile + '.nc'
    else:
        ofile = args[1]
    
    format_options = options.format.split(',')
    file_format = format_options.pop(0)
    f = eval(file_format)(ifile, *format_options)
    pncgen(f, ofile)