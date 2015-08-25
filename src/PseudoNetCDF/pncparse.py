import os
import sys
from warnings import warn
from argparse import ArgumentParser, Action, RawDescriptionHelpFormatter
from cmaqfiles import *
from camxfiles.Memmaps import *
from camxfiles.Readers import irr as irr_read, ipr as ipr_read
from net_balance import mrgaloft, sum_reader, net_reader, ctb_reader
from _getreader import anyfile
from icarttfiles.ffi1001 import ffi1001, ncf2ffi1001
from geoschemfiles import *
from noaafiles import *
from conventions.ioapi import add_cf_from_ioapi, add_cf_from_wrfioapi
from aermodfiles import *
from PseudoNetCDF import PseudoNetCDFFile
from PseudoNetCDF.netcdf import NetCDFFile
from _getreader import getreaderdict
allreaders = getreaderdict()
globals().update(allreaders)
_readernames = [(k.count('.'), k) for k in getreaderdict().keys()]
_readernames.sort()
_readernames = [k for c, k in _readernames]

try:
    from netCDF4 import Dataset as netcdf, MFDataset
except:
    pass

from sci_var import reduce_dim, mesh_dim, slice_dim, getvarpnc, extract, mask_vals, seqpncbo, pncexpr, stack_files, add_attr, convolve_dim, manglenames, removesingleton, merge

class AggCommaString(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        startv = getattr(namespace, self.dest)
        if startv is None:
            startv = []
        setattr(namespace, self.dest, startv + values.split(','))


def getparser(has_ofile, plot_options = False, interactive = False):
    """getparser produces a parser for PseudoNetCDF
    
Parameters
----------
has_ofile : boolean
            Requires the outpath option and processes existence check
plot_options : boolean, optional
               Adds options for plotting including matplotlibrc, spatial
               overlays, and normalization options, default is False
interactive : boolean, optional
              Adds for interactive option, default is False
         
Returns
-------
parser : ArgumentParser with options that are processable with pncparse

    """
    from PseudoNetCDF.register import readers
    parser = ArgumentParser(description = """PseudoNetCDF Argument Parsing


""", formatter_class=RawDescriptionHelpFormatter)
    parser.epilog = """
Detailed Steps
--------------

PseudoNetCDF has many operations and the order often matters. The order is consistent with the order of options in the formatted help. The default order is summarized as:

1. Open with specified reader (-f)
2. Select subset of variables (-v)
2. Add attributes (-a)
4. Apply masks (--mask)
5. Add conventions to support later operations (--to-convention, --from-convention)
6. Combine files via stacking on dimensions (--stack)
7. Slice dimensions (-s --slice)
8. Reduce dimensions (-r --reduce)
9. Convolve dimensions (-c)
10. Extract specific coordinates (--extract)
11. Remove singleton dimensions (--remove-singleton)
12. Apply expressions (--expr then --exprscripts)
13. Apply binary operators (--op_typ)

To impose your own order, use standard options (global options) and then use -- to force positional interpretation of remaining options. In remaining options, use --sep to separate groups of files and options to be evaluated before any global operations."""
    parser.add_argument('ifile', nargs='+', help='path to a file formatted as type -f')
    
    parser.add_argument("-f", "--format", dest = "format", default = 'netcdf', help = "File format (default netcdf); " + '; '.join(_readernames))

    parser.add_argument("--help-format", dest = "helpformat", default = None, help = "Show help for file format (must be one of" + '; '.join(_readernames))
    
    parser.add_argument("--sep", dest = "separator", default = None, help = "Used to separate groups of arguments for parsing (e.g., pncgen -- [options1] file(s)1 [--sep [options2] file(s)2 [... [--sep [optionsN] file(s)N]] ")

    parser.add_argument("--inherit", dest="inherit", action = "store_true", default=False, help = "Allow subparsed sections (separated with -- and --sep) to inherit from global options (-f, --format is always inherited).")

    parser.add_argument("--mangle", dest = "mangle", action = "store_true", default = False, help = "Remove non-standard ascii from names")
    
    parser.add_argument("--remove-singleton", dest = "removesingleton", type = lambda x: [k for k in x.split() if k != ''], default = False, help = "Remove singleton (length 1) dimensions")
    
    # Only has output file if pncgen is called.
    if has_ofile:
        parser.add_argument("-O", "--clobber", dest = "clobber", action = 'store_true', default = False, help = "Overwrite existing file if necessary.")

        parser.add_argument("--out-format", dest = "outformat", default = "NETCDF3_CLASSIC", metavar = "OFMT", help = "File output format (e.g., NETCDF3_CLASSIC, NETCDF4_CLASSIC, NETCDF4;pncgen only)", type = str, choices = 'NETCDF3_CLASSIC NETCDF4_CLASSIC NETCDF4 csv python height_pressure cloud_rain humidity landuse lateral_boundary temperature vertical_diffusivity wind uamiv point_source bpch ffi1001'.split())

        parser.add_argument("--mode", dest = "mode", type = str, default = "w", help = "File mode for writing (w, a or r+ or with unbuffered writes ws, as, or r+s; pncgen only).", choices = 'w a r+ ws as r+s'.split())

        parser.add_argument("--verbose", dest="verbose", action = "store_true", default=False, help = "Provides verbosity with pncgen")

        parser.add_argument('outpath', type = str, help='path to a output file formatted as --out-format')

    else:
        parser.add_argument("-H", "--header", dest="header", action = "store_true", default=False)
        
        parser.add_argument("--full-indices", dest="full_indices",default=None, metavar = "[c|f]", choices = ['c', 'f'])

        parser.add_argument("-l", "--length", dest="line_length", type = int, default=80, metavar = "LEN", help = "CDL line length (pncdump only)")

        parser.add_argument("--float-precision", dest="float_precision", type = int, default=8, metavar = "FDIG", help = "single precision digitis (default 8; pncdump only)")

        parser.add_argument("--double-precision", dest="double_precision", type = int, default=16, metavar = "PDIG", help = "pdig double precision digits (default 16; pncdump only)")

    parser.add_argument("--dump-name", dest = "cdlname", type = str, default = None, help = "Name for display in ncdump")

    parser.add_argument("--coordkeys", dest = "coordkeys", type = str, action = AggCommaString, default = "time time_bounds TFLAG ETFLAG latitude latitude_bounds longitude longitude_bounds lat lat_bnds lon lon_bnds etam_pressure etai_pressure layer_bounds layer47 layer".split(), metavar = "key1,key2", 
                        help = "Variables to be ignored in pncbo.")

    parser.add_argument("-v", "--variables", dest = "variables", default = None, action = AggCommaString,
                        metavar = 'varname1[,varname2[,...,varnameN]',
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")

    parser.add_argument("-a", "--attribute", dest = "attribute", type = str, action = "append", default = [], metavar = "att_nm,var_nm,mode,att_typ,att_val", 
                        help = "Variables have attributes that can be added following nco syntax (--attribute att_nm,var_nm,mode,att_typ,att_val); mode = a,c,d,m,o and att_typ = f,d,l,s,c,b; att_typ is any valid numpy type.")

    parser.add_argument("-m", "--mask", dest = "masks", action = "append", default = [],
                        help = "Masks to apply (e.g., greater,0 or less,0 or values,0, or where,(time[:]%%24<12)[:,None].repeat(10,1))")
        
    parser.add_argument("--from-convention", dest = "fromconv", type = str, default = None, help = "From convention currently only support ioapi")

    parser.add_argument("--to-convention", dest = "toconv", type = str, default = "cf", help = "To convention currently only supports cf")

    parser.add_argument("--stack", dest = "stack", type = str, help = "Concatentate (stack) files on the dimension.")
    
    parser.add_argument("--merge", dest = "merge", action = 'store_true', help = "Combine variables into one file")
    
    parser.add_argument("-s", "--slice", dest = "slice", type = str, action = "append", default = [], metavar = 'dim,start[,stop[,step]]',
                        help = "Variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")

    parser.add_argument("-r", "--reduce", dest = "reduce", type = str, action = "append", default = [], metavar = 'dim,function[,weight]', help = "Variable dimensions can be reduced using dim,function,weight syntax (e.g., --reduce=layer,mean,weight). Weighting is not fully functional.")

    parser.add_argument("--mesh", dest = "mesh", type = str, action = "append", default = [], metavar = 'dim,weight,function', help = "Variable dimensions can be meshed using dim,function,weight syntax (e.g., --mesh=time,0.5,mean).")
    
    parser.add_argument("-c", "--convolve", dest = "convolve", type = str, action = "append", default = [], metavar = 'dim,mode,wgt1,wgt2,...wgtN', help = "Variable dimension is reduced by convolve function (dim,mode,wgt1,wgt2,...wgtN)")    

    parser.add_argument("-e", "--extract", dest = "extract", action = "append", default = [],
                        help = "lon/lat coordinates to extract lon1,lat1/lon2,lat2/lon3,lat3/.../lonN,latN")

    parser.add_argument("--extractmethod", dest = "extractmethod", type = str, default = 'nn', choices = ['nn', 'linear', 'cubic', 'quintic', 'KDTree'],
                        help = "Method for extraction")

    
    parser.add_argument("--op-typ", dest = "operators", type = str, action = 'append', default = [], help = "Operator for binary file operations. Binary file operations use the first two files, then the result and the next file, etc. Use " + " or ".join(['//', '<=', '%%', 'is not', '>>', '&', '==', '!=', '+', '*', '-', '/', '<', '>=', '**', '>', '<<', '|', 'is', '^']))

    parser.add_argument("--expr", dest = "expressions", type = str, action = 'append', default = [], help = "Generic expressions to execute in the context of the file.")

    parser.add_argument("--exprscript", dest = "expressionscripts", type = str, action = 'append', default = [], help = "Generic expressions to execute in the context of the file.")
    if plot_options:
        parser.add_argument("--matplotlibrc", dest = "matplotlibrc", type = str, action = 'append', default = [], help = 'rc options for matplotlib')
        
        parser.add_argument("--plot-commands", dest = "plotcommands", type = str, action = 'append', default = [], help = "Plotting functions to call for all variables expressions to execute in the context of the file.")

        parser.add_argument("--figformat", dest = "figformat", type = str, default = 'png', help = "Any format supported by matplotlib")

        parser.add_argument("--norm", dest = "normalize", type = str, default = None, help = "Typical examples Normalize(), LogNorm(), BoundaryNorm([0, 10, 20, 30, 40], ncolors = 256)")

        parser.add_argument("--colorbar-formatter", dest = "colorbarformatter", type = str, default = None, help = "Typical examples LogFormatter(labelOnlyBase = False), ScalarFormatter(), '%%3g'")

        parser.add_argument("--coastlines", dest = "coastlines", type = bool, default = True, help = "Disable coastlines by setting equal to False")

        parser.add_argument("--countries", dest = "countries", type = bool, default = True, help = "Disable countries by setting equal to False")

        parser.add_argument("--states", dest = "states", type = bool, default = False, help = "Enable states by setting equal to True")

        parser.add_argument("--counties", dest = "counties", type = bool, default = False, help = "Enable counties by setting equal to True")

        parser.add_argument("--shapefiles", dest = "shapefiles", type = str, action = 'append', default = [], help = "Enable custom shapefiles (must be lon, lat)")

        parser.add_argument("--overlay", dest = "overlay", type = bool, default = False, help = "Enable overlay by setting equal to True")

    if interactive:
        parser.add_argument("-i", "--interactive", dest = "interactive", action = 'store_true', default = False, help = "Use interactive mode")
    return parser
        
        
def pncparse(has_ofile = False, plot_options = False, interactive = False, args = None, parser = None):
    """
Parameters
----------
has_ofile : boolean, optional
            Requires the outpath option and processes existence check
            default is False
plot_options : boolean, optional
               Processes matplotlib options before loading matplotlib
               (preprocessing important for backend), default is False.
interactive : boolean, optional
              Only relevant if parser is not provided (see getparser), 
              default is False.
args : list or string, optional
       args are usually taken from the command-line, but can be provided
       in a function call, default is None.
parser : AgumentParser object, optional
         pncparser parser, default getparser(has_ofile, 
                                             plot_options = plot_options, 
                                             interactive = interactive)
         
Returns
-------
(ifiles, args)

ifiles : list of PseudoNetCDFFiles
args : args as parsed

    """
    
    if parser is None:
        parser = getparser(has_ofile, plot_options = plot_options, interactive = interactive)
    helpparser = ArgumentParser(add_help = False)
    helpparser.add_argument("--help", dest = "help", action = 'store_true')
    helpparser.add_argument("--help-format", dest = "helpformat", default = None, help = "Show help for file format (must be one of" + '; '.join(_readernames))
    helpargs, dum = helpparser.parse_known_args(args = args)
    if not helpargs.helpformat is None:
        file_format = helpargs.helpformat.split(',')[0]
        print 'All formats require a "path". Some formats also '
        print 'take extra arguments. All arguments other than '
        print 'the input path must be specified using keyword'
        print 'arguments.'
        print 
        helpformat = allreaders[file_format]
        try:
            import inspect
            print 'Example:'
            idef = inspect.getargspec(helpformat.__init__)
            args = idef.args[2:]
            if idef.defaults is None:
                defs = ['<VAL>'] * len(args)
            else:
                defs = (len(args) - len(idef.defaults)) * ['<VAL>'] + list(idef.defaults)
            longform = 'pncdump -f ' + ','.join([file_format] + [karg + "=%s" % kdef for karg, kdef in zip(args, defs)]) + ' path'
            shortform = 'pncdump -f ' + ','.join([file_format] + [karg + "=%s" % kdef for karg, kdef in zip(args, defs) if kdef == '<VAL>']) + ' path'
            print shortform
            if longform != shortform:
                print ''
                print 'Extended example with keywords:'
                print longform
                print ''
                print '* Any keyword with a non "<VAL>" default can be omitted.'
            print 
        except Exception as e:
            pass
            
        if 'y' == raw_input('Hit Y/y to see detailed help\n').lower():
            help(helpformat)
        exit()
    subparser = getparser(has_ofile = False, plot_options = plot_options, interactive = interactive)
    args = parser.parse_args(args = args)
    subargs = split_positionals(subparser, args)
    ifiles = []
    ipaths = []
    for subarg in subargs:
        ipaths.extend(subarg.ifile)
        ifiles.extend(pncprep(subarg)[0])
    args.ifile = ifiles
    args.ipath = ipaths
    #if args.stack is not None:
    #    pass
    #elif len(args.ifile) == 0 or (len(args.ifile) - len(args.operators) - has_ofile) > 1:
    #    print 'Too many arguments', len(args.ifile), len(args.operators), has_ofile
    #    parser.print_help()
    #    exit()
    if has_ofile:
        if not args.clobber and os.path.exists(args.outpath):
            parser.error(message = 'Output path (%s) exists; enable clobber (--clobber or -O) to overwrite.' % (args.outpath,))
    else:
        args.outpath = None

    if plot_options:
        from matplotlib import rcParams
        for rcassign in args.matplotlibrc:
            key, value = rcassign.split('=')
            rcParams[key] = value

    return pncprep(args)

def split_positionals(parser, args):
    positionals = args.ifile
    parser.set_defaults(**dict([(k, v) for k, v in args._get_kwargs() if args.inherit or k == 'format']))
    outs = []
    last_split = 0
    for i in range(len(positionals)):
        if positionals[i] in ('--', '--sep'):
            outs.append(positionals[last_split:i])
            last_split = i + 1
    outs.append(positionals[last_split:])        
    outs = map(parser.parse_args, outs)
    
    for out in outs:
        if out.cdlname is None:
            out.cdlname = ', '.join(out.ifile)
    if len(outs) == 1 and args.cdlname is None:
        args.cdlname = outs[0].cdlname
    return outs


def pncprep(args):
    #nifiles = len(args.ifile) - has_ofile
    ipaths = args.ifile[:]
    fs = getfiles(ipaths, args)
    fs = subsetfiles(fs, args)
    fs = seqpncbo(args.operators, fs, coordkeys = args.coordkeys)
    for expr in args.expressions:
        fs = [pncexpr(expr, f) for f in fs]
    for script in args.expressionscripts:
        expr = open(script).read()
        fs = [pncexpr(expr, f) for f in fs]
    return fs, args

def subsetfiles(ifiles, args):
    fs = []
    for f in ifiles:
        for opts in args.slice:
            f = slice_dim(f, opts)
        for opts in args.reduce:
            f = reduce_dim(f, opts, metakeys = args.coordkeys)
        for opts in args.mesh:
            f = mesh_dim(f, opts)
        for opts in args.convolve:
            f = convolve_dim(f, opts)
        if len(args.extract) > 0:
            f = extract(f, args.extract, method = args.extractmethod)
        if args.removesingleton != False:
            f = removesingleton(f)
        fs.append(f)
    return fs

def getfiles(ipaths, args):
    fs = []
    for ipath in ipaths:
        format_options = args.format.split(',')
        file_format = format_options.pop(0)
        format_options = eval('dict(' + ', '.join(format_options) + ')')
        if isinstance(ipath, (PseudoNetCDFFile, NetCDFFile)):
            f = ipath
        elif isinstance(ipath, (str, unicode)) :
            try:
                if file_format in allreaders:
                    f = allreaders[file_format](ipath, **format_options)
                else:
                    f = eval(file_format)(ipath, **format_options)
            except Exception as e:
                raise IOError('Unable to open path with %s(path, **%s)\n\tpath="%s"\n\terror="%s"' % (file_format, str(format_options), ipath, str(e)))
        else:
            warn('File is type %s, which is unknown' % type(ipath))
            f = ipath
        history = getattr(f, 'history', '')
        history += ' '.join(sys.argv) + ';'
        laddconv = args.fromconv is not None and args.toconv is not None
        lslice = len(args.slice + args.reduce) > 0
        lexpr = len(args.expressions) > 0
        if args.variables is not None:
            f = getvarpnc(f, args.variables, coordkeys = args.coordkeys)
        elif laddconv or lslice or lexpr:
            f = getvarpnc(f, None)
        for opts in args.attribute:
            add_attr(f, opts)
        for opts in args.masks:
            f = mask_vals(f, opts, metakeys = args.coordkeys)
        if laddconv:
            try:
                eval('add_%s_from_%s' % (args.toconv, args.fromconv))(f)
            except Exception as e:
                warn('Cannot add %s from %s; %s' % (args.toconv, args.fromconv, str(e)))
        if args.cdlname is None:
            args.cdlname = ipath
        
        try:
            setattr(f, 'history', history)
        except: pass
        if args.mangle:
            f = manglenames(f)
        fs.append(f)
    if args.stack is not None:
        fs = [stack_files(fs, args.stack, coordkeys = args.coordkeys)]
    if args.merge:
        fs = [merge(fs)]
    return fs
