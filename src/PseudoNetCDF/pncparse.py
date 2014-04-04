import os
import sys
from _getreader import getreaderdict
globals().update(getreaderdict())
from argparse import ArgumentParser, Action
from cmaqfiles import *
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

class AggCommaString(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        startv = getattr(namespace, self.dest)
        if startv is None:
            startv = []
        setattr(namespace, self.dest, startv + values.split(','))


def pncparser(has_ofile):
    parser = ArgumentParser(description = 'PseudoNetCDF Argument Parsing')
    parser.add_argument('ifile', nargs='+', help='path to a file formatted as type -f')
    parser.add_argument("-f", "--format", dest = "format", default = 'netcdf', help = "File format (default netcdf)")

    # Only has output file if pncgen is called.
    if has_ofile:
        parser.add_argument('outpath', type = str, help='path to a output file formatted as --out-format')
        parser.add_argument("-O", "--clobber", dest = "clobber", action = 'store_true', default = False, help = "Overwrite existing file if necessary.")
        parser.add_argument("--out-format", dest = "outformat", default = "NETCDF3_CLASSIC", metavar = "OFMT", help = "File output format (e.g., NETCDF3_CLASSIC, NETCDF4_CLASSIC, NETCDF4;pncgen only)", type = str, choices = 'NETCDF3_CLASSIC NETCDF4_CLASSIC NETCDF4 height_pressure cloud_rain humidity landuse lateral_boundary temperature vertical_diffusivity wind uamiv point_source bpch'.split())
        parser.add_argument("--verbose", dest="verbose", action = "store_true", default=False, help = "Provides verbosity with pncgen")

        parser.add_argument("--mode", dest = "mode", type = str, default = "w", help = "File mode for writing (w, a or r+ or with unbuffered writes ws, as, or r+s; pncgen only).", choices = 'w a r+ ws as r+s'.split())

    else:
        parser.add_argument("-H", "--header", dest="header", action = "store_true", default=False)
        
        parser.add_argument("--full-indices", dest="full_indices",default=None, metavar = "[c|f]", choices = ['c', 'f'])

        parser.add_argument("-l", "--length", dest="line_length", type = int, default=80, metavar = "LEN", help = "CDL line length (pncdump only)")

        parser.add_argument("--float-precision", dest="float_precision", type = int, default=8, metavar = "FDIG", help = "single precision digitis (default 8; pncdump only)")

        parser.add_argument("--double-precision", dest="double_precision", type = int, default=16, metavar = "PDIG", help = "pdig double precision digits (default 16; pncdump only)")

    parser.add_argument("--dump-name", dest = "cdlname", type = str, default = None, help = "Name for display in ncdump")

    parser.add_argument("-v", "--variables", dest = "variables", default = None, action = AggCommaString,
                        metavar = 'varname1[,varname2[,...,varnameN]',
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")

    parser.add_argument("--from-convention", dest = "fromconv", type = str, default = None, help = "From convention currently only support ioapi")

    parser.add_argument("--to-convention", dest = "toconv", type = str, default = None, help = "To convention currently only supports cf")

    parser.add_argument("-s", "--slice", dest = "slice", type = str, action = "append", default = [], metavar = 'dim,start[,stop[,step]]',
                        help = "Variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")

    parser.add_argument("-r", "--reduce", dest = "reduce", type = str, action = "append", default = [], metavar = 'dim,function[,weight]', help = "Variable dimensions can be reduced using dim,function,weight syntax (e.g., --reduce=layer,mean,weight). Weighting is not fully functional.")

    parser.add_argument("-c", "--convolve", dest = "convolve", type = str, action = "append", default = [], metavar = 'dim,mode,wgt1,wgt2,...wgtN', help = "Variable dimension is reduced by convolve function (dim,mode,wgt1,wgt2,...wgtN)")
    
    parser.add_argument("-a", "--attribute", dest = "attribute", type = str, action = "append", default = [], metavar = "att_nm,var_nm,mode,att_typ,att_val", 
                        help = "Variables have attributes that can be added following nco syntax (--attribute att_nm,var_nm,mode,att_typ,att_val); mode = a,c,d,m,o and att_typ = f,d,l,s,c,b; att_typ is any valid numpy type.")

    parser.add_argument("--post-attribute", dest = "postattribute", type = str, action = "append", default = [], metavar = "att_nm,var_nm,mode,att_typ,att_val", 
                        help = "Variables have attributes that can be added following nco syntax (--attribute att_nm,var_nm,mode,att_typ,att_val); mode = a,c,d,m,o and att_typ = f,d,l,s,c,b; att_typ is any valid numpy type.")

    parser.add_argument("--coordkeys", dest = "coordkeys", type = lambda c_: str(c_).split(), action = "append", default = ["time time_bounds latitude latitude_bounds longitude longitude_bounds lat lat_bnds lon lon_bnds".split()], metavar = "key1,key2", 
                        help = "Variables to be ignored in pncbo.")


    parser.add_argument("--mesh", dest = "mesh", type = str, action = "append", default = [], metavar = 'dim,weight,function', help = "Variable dimensions can be reduced using dim,function,weight syntax (e.g., --mesh=time,0.5,mean).")
    

    parser.add_argument("-e", "--extract", dest = "extract", action = "append", default = [],
                        help = "lon/lat coordinates to extract lon1,lat1/lon2,lat2/lon3,lat3/.../lonN,latN")

    parser.add_argument("-m", "--mask", dest = "masks", action = "append", default = [],
                        help = "Masks to apply (e.g., greater,0 or less,0 or values,0)")
        
    
    parser.add_argument("--op-typ", dest = "operators", type = str, action = 'append', default = [], help = "Operator for binary file operations. Binary file operations use the first two files, then the result and the next file, etc. Use " + " or ".join(['//', '<=', '%%', 'is not', '>>', '&', '==', '!=', '+', '*', '-', '/', '<', '>=', '**', '>', '<<', '|', 'is', '^']))

    parser.add_argument("--stack", dest = "stack", type = str, help = "Concatentate (stack) files on the dimension.")
    
    parser.add_argument("--op-first", dest = "operatorsfirst", action = 'store_true', default = False, help = "Use operations before slice/aggregate.")

    parser.add_argument("--expr", dest = "expressions", type = str, action = 'append', default = [], help = "Generic expressions to execute in the context of the file.")


    args = parser.parse_args()
    args.coordkeys = reduce(list.__add__, args.coordkeys)
    #if args.stack is not None:
    #    pass
    #elif len(args.ifile) == 0 or (len(args.ifile) - len(args.operators) - has_ofile) > 1:
    #    print 'Too many arguments', len(args.ifile), len(args.operators), has_ofile
    #    parser.print_help()
    #    exit()
    
    #nifiles = len(args.ifile) - has_ofile
    ipaths = args.ifile[:]
    if has_ofile:
        if not args.clobber and os.path.exists(args.outpath):
            parser.error(message = 'Output path (%s) exists; enable clobber (--clobber or -O) to overwrite.' % (args.outpath,))
    else:
        args.outpath = None
    
    fs = getfiles(ipaths, args)
    
    if args.operatorsfirst: fs = seqpncbo(args.operators, fs, coordkeys = args.coordkeys)
    fs = subsetfiles(fs, args)
    if not args.operatorsfirst: fs = seqpncbo(args.operators, fs, coordkeys = args.coordkeys)
    for expr in args.expressions:
        for f in fs:
            f = pncexpr(expr, f)
    decorate(fs, args)
    return fs, args

def decorate(ifiles, args):
    for f in ifiles:
        for opts in args.postattribute:
            add_attr(f, opts)

def subsetfiles(ifiles, args):
    fs = []
    for f in ifiles:
        for opts in args.slice:
            f = slice_dim(f, opts)
        for opts in args.reduce:
            f = reduce_dim(f, opts)
        for opts in args.mesh:
            f = mesh_dim(f, opts)
        for opts in args.convolve:
            f = convolve_dim(f, opts)
        if len(args.extract) > 0:
            f = extract(f, args.extract)
        fs.append(f)
    return fs

def getfiles(ipaths, args):
    fs = []
    for ipath in ipaths:
        format_options = args.format.split(',')
        file_format = format_options.pop(0)
        format_options = eval('dict(' + ', '.join(format_options) + ')')
        f = eval(file_format)(ipath, **format_options)
        history = getattr(f, 'history', '')
        history += ' '.join(sys.argv) + ';'
        laddconv = args.fromconv is not None and args.toconv is not None
        lslice = len(args.slice + args.reduce) > 0
        lexpr = len(args.expressions) > 0
        if args.variables is not None:
            f = getvarpnc(f, args.variables)
        elif laddconv or lslice or lexpr:
            f = getvarpnc(f, None)
        for opts in args.attribute:
            add_attr(f, opts)
        for opts in args.masks:
            f = mask_vals(f, opts)
        if laddconv:
            eval('add_%s_from_%s' % (args.toconv, args.fromconv))(f)
        if args.cdlname is None:
            args.cdlname = ipath
        try:
            setattr(f, 'history', history)
        except: pass
        fs.append(f)
    if args.stack is not None:
        fs = [stack_files(fs, args.stack)]
    return fs
