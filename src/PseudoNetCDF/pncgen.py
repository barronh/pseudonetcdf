from types import MethodType
from sci_var import Pseudo2NetCDF
from netcdf import NetCDFFile

def pncgen(f,outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4'):
    p2n = Pseudo2NetCDF()
    return p2n.convert(f, outpath, inmode = inmode, outmode = outmode, format = format)
    
if __name__ == '__main__':
    from optparse import OptionParser
    from camxfiles.Memmaps import *
    from icarttfiles.ffi1001 import ffi1001
    from geoschemfiles import bpch
    from sci_var import reduce_dim, slice_dim, getvarpnc

    parser = OptionParser()
    parser.set_usage("""Usage: python -m %prog [-f uamiv|bpch|ffi1001|...] ifile [ofile]

    ifile - path to a file formatted as type -f
    ofile - path to the desired output
    
    """)

    parser.add_option("-f", "--format", dest = "format", default = 'uamiv', help = "File format")
    
    parser.add_option("-v", "--variables", dest = "variables", default = None,
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")

    parser.add_option("-s", "--slice", dest = "slice", type = "string", action = "append", default = [],
                        help = "bpch variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")
    
    parser.add_option("-r", "--reduce", dest = "reduce", type = "string", action = "append", default = [], help = "bpch variable dimensions can be reduced using dim,function,weight syntax (e.g., --reduce=layer,mean,weight). Weighting is not fully functional.")


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
    format_options = eval('dict(' + ', '.join(format_options) + ')')
    f = eval(file_format)(ifile, **format_options)
    
    if options.variables is not None:
        f = getvarpnc(f, options.variables.split(','))
    for opts in options.slice:
        f = slice_dim(f, opts)
    for opts in options.reduce:
        f = reduce_dim(f, opts)

    pncgen(f, ofile)