import sys
from _getreader import getreaderdict
globals().update(getreaderdict())
from optparse import OptionParser
from camxfiles.Memmaps import *
from camxfiles.Readers import irr as irr_read, ipr as ipr_read
from net_balance import mrgaloft, sum_reader, net_reader, ctb_reader
from _getreader import anyfile
from icarttfiles.ffi1001 import ffi1001
from geoschemfiles import *
from noaafiles import *
from conventions.ioapi import add_cf_from_ioapi

try:
    from netCDF4 import Dataset as netcdf
except:
    pass
from sci_var import reduce_dim, mesh_dim, slice_dim, getvarpnc, extract, mask_vals, seqpncbo, pncexpr, stack_files, add_attr, convolve_dim

def pncparser(has_ofile):
    parser = OptionParser()
    parser.set_usage("""Usage: python -m %prog [-f netcdf|uamiv|bpch|ffi1001|...] ifile [ifile2 [ifile3 [...]]] [ofile]

    ifile - path to a file formatted as type -f
    ofile - path to the desired output (pncgen only)
    
    """)

    parser.add_option("-f", "--format", dest = "format", default = 'netcdf', help = "File format (default netcdf)")

    parser.add_option("-H", "--header", dest="header", action = "store_true", default=False)

    parser.add_option("-v", "--variables", dest = "variables", default = None,
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")

    parser.add_option("", "--from-convention", dest = "fromconv", type = "string", default = None, help = "From convention currently only support ioapi")

    parser.add_option("", "--to-convention", dest = "toconv", type = "string", default = None, help = "To convention currently only supports cf")

    parser.add_option("-s", "--slice", dest = "slice", type = "string", action = "append", default = [],
                        help = "Variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")

    parser.add_option("-r", "--reduce", dest = "reduce", type = "string", action = "append", default = [], help = "Variable dimensions can be reduced using dim,function,weight syntax (e.g., --reduce=layer,mean,weight). Weighting is not fully functional.")

    parser.add_option("-c", "--convolve", dest = "convolve", type = "string", action = "append", default = [], help = "Variable dimension is reduced by convolve function (dim,mode,wgt1,wgt2,...wgtN)")
    
    parser.add_option("-a", "--attribute", dest = "attribute", type = "string", action = "append", default = [],
                        help = "Variables have attributes that can be added following nco syntax (--attribute att_nm,var_nm,mode,att_typ,att_val); mode = a,c,d,m,o and att_typ = f,d,l,s,c,b; att_typ is any valid numpy type.")

    parser.add_option("", "--post-attribute", dest = "postattribute", type = "string", action = "append", default = [],
                        help = "Variables have attributes that can be added following nco syntax (--attribute att_nm,var_nm,mode,att_typ,att_val); mode = a,c,d,m,o and att_typ = f,d,l,s,c,b; att_typ is any valid numpy type.")


    parser.add_option("", "--mesh", dest = "mesh", type = "string", action = "append", default = [], help = "Variable dimensions can be reduced using dim,function,weight syntax (e.g., --mesh=time,0.5,mean).")
    

    parser.add_option("-e", "--extract", dest = "extract", action = "append", default = [],
                        help = "lon/lat coordinates to extract")

    parser.add_option("-m", "--mask", dest = "masks", action = "append", default = [],
                        help = "Masks to apply (e.g., greater,0 or less,0 or values,0)")
    
    parser.add_option("", "--full-indices", dest="full_indices",default=None, metavar = "[c|f]", choices = ['c', 'f'])

    parser.add_option("-l", "--length", dest="line_length", type = "int", default=80, metavar = "LEN", help = "CDL line length (pncdump only)")

    parser.add_option("", "--verbose", dest="verbose", action = "store_true", default=False, help = "Provides verbosity with pncgen")

    parser.add_option("", "--float-precision", dest="float_precision", type="int", default=8, metavar = "FDIG", help = "single precision digitis (default 8; pncdump only)")

    parser.add_option("", "--double-precision", dest="double_precision", type="int", default=16, metavar = "PDIG", help = "pdig double precision digits (default 16; pncdump only)")

    parser.add_option("", "--dump-name", dest = "cdlname", type = "string", default = None, help = "Name for display in ncdump")

    parser.add_option("", "--op-typ", dest = "operators", type = "string", action = 'append', default = [], help = "Operator for binary file operations. Binary file operations use the first two files, then the result and the next file, etc. Use " + " or ".join(['//', '<=', '%', 'is not', '>>', '&', '==', '!=', '+', '*', '-', '/', '<', '>=', '**', '>', '<<', '|', 'is', '^']))

    parser.add_option("", "--stack", dest = "stack", type = "string", help = "Concatentate (stack) files on the dimension.")
    
    parser.add_option("", "--op-first", dest = "operatorsfirst", action = 'store_true', default = False, help = "Use operations before slice/aggregate.")

    parser.add_option("", "--expr", dest = "expressions", type = "string", action = 'append', default = [], help = "Generic expressions to execute in the context of the file.")
    
    parser.add_option("", "--out-format", dest = "outformat", default = "NETCDF3_CLASSIC", metavar = "OFMT", help = "File output format (e.g., NETCDF3_CLASSIC, NETCDF4_CLASSIC, NETCDF4;pncgen only)", type = "choice", choices = 'NETCDF3_CLASSIC NETCDF4_CLASSIC NETCDF4 height_pressure cloud_rain humidity landuse lateral_boundary temperature vertical_diffusivity wind uamiv point_source bpch'.split())

    parser.add_option("", "--mode", dest = "mode", type = "choice", default = "w", help = "File mode for writing (w, a or r+ or with unbuffered writes ws, as, or r+s; pncgen only).", choices = 'w a r+ ws as r+s'.split())

    (options, args) = parser.parse_args()
    
    if hasattr(options, 'stack'):
        pass
    elif len(args) == 0 or (len(args) - len(options.operators) - has_ofile) > 1:
        print 'Too many arguments', len(args), len(options.operators), has_ofile
        parser.print_help()
        exit()
    
    nifiles = len(args) - has_ofile
    ipaths = args[:nifiles]
    if has_ofile:
        if len(args) == 1:
            ofile = ifiles[0] + '.nc'
        else:
            ofile = args[-1]
    else:
        ofile = None
    fs = getfiles(ipaths, options)
    
    if options.operatorsfirst: fs = seqpncbo(options.operators, fs)
    fs = subsetfiles(fs, options)
    if not options.operatorsfirst: fs = seqpncbo(options.operators, fs)
    f, = fs
    for expr in options.expressions:
        f = pncexpr(expr, f)
    decorate(fs, options)
    return f, ofile, options

def decorate(ifiles, options):
    for f in ifiles:
        for opts in options.postattribute:
            add_attr(f, opts)

def subsetfiles(ifiles, options):
    fs = []
    for f in ifiles:
        for opts in options.slice:
            f = slice_dim(f, opts)
        for opts in options.reduce:
            f = reduce_dim(f, opts)
        for opts in options.mesh:
            f = mesh_dim(f, opts)
        for opts in options.convolve:
            f = convolve_dim(f, opts)
        if len(options.extract) > 0:
            extract(f, options.extract)
        fs.append(f)
    return fs

def getfiles(ipaths, options):
    fs = []
    for ipath in ipaths:
        format_options = options.format.split(',')
        file_format = format_options.pop(0)
        format_options = eval('dict(' + ', '.join(format_options) + ')')
        f = eval(file_format)(ipath, **format_options)
        history = getattr(f, 'history', '')
        history += ' '.join(sys.argv) + ';'
        laddconv = options.fromconv is not None and options.toconv is not None
        lslice = len(options.slice + options.reduce) > 0
        if options.variables is not None:
            f = getvarpnc(f, options.variables.split(','))
        elif laddconv or lslice:
            f = getvarpnc(f, None)
        for opts in options.attribute:
            add_attr(f, opts)
        for opts in options.masks:
            f = mask_vals(f, opts)
        if laddconv:
            eval('add_%s_from_%s' % (options.toconv, options.fromconv))(f)
        if options.cdlname is None:
            options.cdlname = ipath
        try:
            setattr(f, 'history', history)
        except: pass
        fs.append(f)
    if options.stack is not None:
        fs = stack_files(fs, options.stack)
    return fs
