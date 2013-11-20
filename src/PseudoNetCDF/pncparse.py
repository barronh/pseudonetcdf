import sys
from optparse import OptionParser
from camxfiles.Memmaps import *
from camxfiles.Readers import irr as irr_read, ipr as ipr_read
from icarttfiles.ffi1001 import ffi1001
from geoschemfiles import *
from conventions.ioapi import add_cf_from_ioapi
try:
    from netCDF4 import Dataset as netcdf
except:
    pass
from sci_var import reduce_dim, slice_dim, getvarpnc, extract, mask_vals

def pncparser():
    parser = OptionParser()
    parser.set_usage("""Usage: python -m %prog [-f netcdf|uamiv|bpch|ffi1001|...] ifile [ofile]

    ifile - path to a file formatted as type -f
    ofile - path to the desired output (pncgen only)
    
    """)

    parser.add_option("-f", "--format", dest = "format", default = 'netcdf', help = "File format (default netcdf)")

    parser.add_option("-H", "--header", dest="header", action = "store_true", default=False)

    parser.add_option("", "--variables", dest = "variables", default = None,
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")

    parser.add_option("", "--from-convention", dest = "fromconv", type = "string", default = None, help = "From convention currently only support ioapi")

    parser.add_option("", "--to-convention", dest = "toconv", type = "string", default = None, help = "To convention currently only supports cf")

    parser.add_option("-s", "--slice", dest = "slice", type = "string", action = "append", default = [],
                        help = "bpch variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")

    parser.add_option("-r", "--reduce", dest = "reduce", type = "string", action = "append", default = [], help = "bpch variable dimensions can be reduced using dim,function,weight syntax (e.g., --reduce=layer,mean,weight). Weighting is not fully functional.")

    parser.add_option("-e", "--extract", dest = "extract", action = "append", default = [],
                        help = "lon/lat coordinates to extract")

    parser.add_option("-m", "--mask", dest = "masks", action = "append", default = [],
                        help = "Masks to apply (e.g., greater,0 or less,0 or values,0)")
    
    parser.add_option("", "--full-indices", dest="full_indices",default=None, metavar = "[C|F]")

    parser.add_option("-l", "--length", dest="line_length", type = "int", default=80, metavar = "LEN", help = "CDL line length (pncdump only)")

    parser.add_option("-v", "--verbose", dest="verbose", action = "store_true", default=False, help = "Provides verbosity with pncgen")

    parser.add_option("", "--float-precision", dest="float_precision", type="int", default=8, metavar = "FDIG", help = "single precision digitis (default 8; pncdump only)")

    parser.add_option("", "--double-precision", dest="double_precision", type="int", default=16, metavar = "PDIG", help = "pdig double precision digits (default 16; pncdump only)")

    parser.add_option("", "--dump-name", dest = "cdlname", type = "string", default = None, help = "Name for display in ncdump")
    
    parser.add_option("", "--out-format", dest = "outformat", default = "NETCDF3_CLASSIC", metavar = "OFMT", help = "File output format (e.g., NETCDF3_CLASSIC, NETCDF4_CLASSIC, NETCDF4;pncgen only)")

    parser.add_option("", "--mode", dest = "mode", type = "string", default = "w", help = "File mode for writing (w, a or r+ or with unbuffered writes ws, as, or r+s; pncgen only).")


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
    history = getattr(f, 'history', '')
    history += ' '.join(sys.argv) + ';'
    laddconv = options.fromconv is not None and options.toconv is not None
    lslice = len(options.slice + options.reduce) > 0
    if options.variables is not None:
        f = getvarpnc(f, options.variables.split(','))
    elif laddconv or lslice:
        f = getvarpnc(f, None)
    for opts in options.masks:
        f = mask_vals(f, opts)
    for opts in options.slice:
        f = slice_dim(f, opts)
    for opts in options.reduce:
        f = reduce_dim(f, opts)
    if laddconv:
        eval('add_%s_from_%s' % (options.toconv, options.fromconv))(f)
    if len(options.extract) > 0:
        extract(f, options.extract)
    if options.cdlname is None:
        options.cdlname = ifile
    try:
        setattr(f, 'history', history)
    except: pass
    return f, ofile, options