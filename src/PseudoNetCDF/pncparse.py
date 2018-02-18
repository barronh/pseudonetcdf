from __future__ import print_function, unicode_literals
import os
import sys
from warnings import warn
from argparse import ArgumentParser, Action, RawDescriptionHelpFormatter
from .cmaqfiles import *
from .camxfiles.Memmaps import *
from .camxfiles.Readers import irr as irr_read, ipr as ipr_read
from .net_balance import mrgaloft, sum_reader, net_reader, ctb_reader
from ._getreader import anyfile, getreaderdict
from ._getwriter import getwriterdict
from .icarttfiles.ffi1001 import ffi1001, ncf2ffi1001
from .geoschemfiles import *
from .noaafiles import *
from .conventions.ioapi import *
from .aermodfiles import *
from PseudoNetCDF import PseudoNetCDFFile, netcdf
from PseudoNetCDF.netcdf import NetCDFFile

allreaders = getreaderdict()
allwriters = getwriterdict()
globals().update(allreaders)
_readernames = [(k.count('.') if k[:1] != '_' else 9999, k) for k in getreaderdict().keys()]
_readernames.sort()
_readernames = [k for c, k in _readernames]
_writernames = [(k.count('.'), k) for k in getwriterdict().keys()]
_writernames.sort()
_writernames = [k for c, k in _writernames]

try:
    from netCDF4 import MFDataset
except:
    pass

from .sci_var import reduce_dim, mesh_dim, slice_dim, getvarpnc, extract, mask_vals, seqpncbo, pncexpr, stack_files, add_attr, convolve_dim, manglenames, removesingleton, merge, extract_from_file, pncrename, WrapPNC

        
from argparse import SUPPRESS

class PNCArgumentParser(ArgumentParser):
    def exit(self, status=0, message=None):
        import sys as _sys
        if message:
            self._print_message(message, _sys.stderr)

class _HelpAction(Action):
    def __init__(self,
                 option_strings,
                 dest=SUPPRESS,
                 default=SUPPRESS,
                 help=None):
        super(_HelpAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        setattr(namespace, self.dest, True)

class _HelpListFormats(Action):
    def __init__(self,
                 option_strings,
                 dest=SUPPRESS,
                 default=SUPPRESS,
                 help=None):
        super(_HelpListFormats, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print('The formats listed below are available for the following options')
        print('-f FORMAT or --format FORMAT or --help-format FORMAT')
        print('where FORMAT is one of the options below')
        print('\t' + '\n\t'.join(_readernames))
        print('**Some readers are listed twice (e.g., without dotted form)')
        setattr(namespace, 'help', True)

class _HelpFormat(Action):
    def __init__(self,
                 option_strings,
                 dest=SUPPRESS,
                 default=SUPPRESS,
                 help=None):
        super(_HelpFormat, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=1,
            help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        file_format = values[0].split(',')[0]
        print('')
        print('All formats require a "path". Some formats also ')
        print('take extra arguments. All arguments other than ')
        print('the input path must be specified using keyword')
        print('arguments.')
        print('')
        helpformat = allreaders[file_format]
        try:
            import inspect
            print('Example:')
            idef = inspect.getargspec(helpformat.__init__)
            args = idef.args[2:]
            if idef.defaults is None:
                defs = ['<VAL>'] * len(args)
            else:
                defs = (len(args) - len(idef.defaults)) * ['<VAL>'] + list(idef.defaults)
            longform = 'pncdump -f ' + ','.join([file_format] + [karg + "=%s" % kdef for karg, kdef in zip(args, defs)]) + ' path'
            shortform = 'pncdump -f ' + ','.join([file_format] + [karg + "=%s" % kdef for karg, kdef in zip(args, defs) if kdef == '<VAL>']) + ' path'
            print(shortform)
            if longform != shortform:
                print('')
                print('Extended example with keywords:')
                print(longform)
                print('')
                print('* Any keyword with a non "<VAL>" default can be omitted.')
            print('')
        except Exception as e:
            print('** Short help not available for ' + file_format)
            print('')
        
        if 'y' == input('Hit Y/y to see detailed help\n').lower():
            help(helpformat)
        setattr(namespace, self.dest, values)
        setattr(namespace, 'help', True)

class AggCommaString(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

_plotcmds = ['plot', 'plot2d', 'plotts', 'plotprofile', 'plotscatter']
_allcmds = ['gen', 'dump', 'eval', 'map'] + _plotcmds

def add_basic_options(parser):
    parser.add_argument('ifiles', nargs='*', help='path to a file formatted as type -f')
    parser.add_argument("--verbose", dest="verbose", action = "count", default=0, help = "Provides verbosity with pncgen")

    try:
        parser.add_argument("--help", dest="help", action = _HelpAction, default=False, help = "Displays help")
    except:
        pass

    parser.add_argument('--pnc', action = 'append', default = [], help='Set of pseudonetcdf commands to be process separately')

    parser.add_argument("-f", "--format", dest = "format", default = 'netcdf', metavar = '{see --list-formats for choices}', help = "File format (default netcdf), can be one of the choices listed, or an expression that evaluates to a reader. Keyword arguments are passed via ,kwd=value.")

    parser.add_argument("--list-formats", dest = "help", default = None, action = _HelpListFormats, help = "Show format options for -f")

    parser.add_argument("--help-format", dest = 'helpformat', action = _HelpFormat, default = None, help = "Show help for file format (must be one of the options for -f)")
    
    parser.add_argument("--sep", dest = "separator", action = 'store_true', default = False, help = "Used to separate groups of arguments for parsing (e.g., pncgen -- [options1] file(s)1 [--sep [options2] file(s)2 [... [--sep [optionsN] file(s)N]] ")

    parser.add_argument("--inherit", dest="inherit", action = "store_true", default=False, help = "Allow subparsed sections (separated with -- and --sep) to inherit from global options (-f, --format is always inherited).")

    parser.add_argument("--mangle", dest = "mangle", action = "store_true", default = False, help = "Remove non-standard ascii from names")

    parser.add_argument("--rename", dest = "rename", action = "append", default = [], help = "Provide pairs of strings to be substituted --rename=type,oldkey,newkey (type: v = variable; d = dimension;)")
    
    parser.add_argument("--remove-singleton", dest = "removesingleton", type = lambda x: [k for k in x.split(',') if k != ''], default = [], help = "Remove singleton (length 1) dimensions")
    
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

    parser.add_argument("--extract-file", dest = "extractfile", action = "append", default = [],
                        help = "pncparse options for file")

    parser.add_argument("--extractmethod", dest = "extractmethod", type = str, default = 'nn', choices = ['ll2ij', 'nn', 'linear', 'cubic', 'quintic', 'KDTree'],
                        help = "Method for extraction")

    
    parser.add_argument("--op-typ", dest = "operators", type = str, action = 'append', default = [], help = "Operator for binary file operations. Binary file operations use the first two files, then the result and the next file, etc. Use " + " or ".join(['//', '<=', '%%', 'is not', '>>', '&', '==', '!=', '+', '*', '-', '/', '<', '>=', '**', '>', '<<', '|', 'is', '^']))

    parser.add_argument("--expr", dest = "expressions", type = str, action = 'append', default = [], help = "Generic expressions to execute in the context of the file.")

    parser.add_argument("--exprscript", dest = "expressionscripts", type = str, action = 'append', default = [], help = "Generic expressions to execute in the context of the file.")

def add_2d_options(parser):
    parser.add_argument('--swapaxes', action = 'store_true', help = 'Swap x-y axes')

def add_eval_options(parser):
    from .pnceval import __all__ as evalfuncs
    add_interactive_options(parser)
    parser.add_argument('--funcs', default = evalfuncs, type = lambda x: x.split(','), help='Functions to evaluate split by , (default: %s)' % ','.join(evalfuncs))

def add_map_options(parser):
    parser.add_argument("--iter", dest = "iter", action = 'append', default = [], help = "Create plots for each element of specified dimension (e.g., --iter=time).")

    parser.add_argument("--resolution", dest = "resolution", choices = ['c', 'i', 'h'], default = 'c', help = "Use coarse (c), intermediate (i), or high (h) resolution for maps.")
    
    parser.add_argument("--no-coastlines", dest = "coastlines", action = 'store_false', help = "Disable coastlines by setting equal to False")

    parser.add_argument("--no-countries", dest = "countries", action = 'store_false', help = "Disable countries by setting equal to False")

    parser.add_argument("--states", dest = "states", action = 'store_true', help = "Enable states by setting equal to True")

    parser.add_argument("--counties", dest = "counties", action = 'store_true', help = "Enable counties by setting equal to True")

    parser.add_argument("--shapefile", "--shapefiles", dest = "shapefiles", type = str, action = 'append', default = [], help = "Enable custom shapefiles (must be lon, lat). Keyword arguments for display can be added with commas (see -f for syntax).")
    
def add_dump_options(parser):
    parser.add_argument("-H", "--header", dest="header", action = "store_true", default=False)
        
    parser.add_argument("-t", "--timestring", dest="timestring", action = "store_true", default=False)
    
    parser.add_argument("--full-indices", dest="full_indices",default=None, metavar = "[c|f]", choices = ['c', 'f'], help = "Provide indices in CDL using either C or Fortran style indexes. C style is 0-based and ordered from slowest iterating dimension to fastest. Fortran style is 1-based and ordered from fastest to slowest iterating dimension")

    parser.add_argument("-l", "--length", dest="line_length", type = int, default=80, metavar = "LEN", help = "CDL line length (pncdump only)")

    parser.add_argument("--float-precision", dest="float_precision", type = int, default=8, metavar = "FDIG", help = "single precision digitis (default 8; pncdump only)")

    parser.add_argument("--double-precision", dest="double_precision", type = int, default=16, metavar = "PDIG", help = "pdig double precision digits (default 16; pncdump only)")

    parser.add_argument("--dump-name", dest = "cdlname", type = str, default = None, help = "Name for display in ncdump")

def add_output_options(parser):
    parser.add_argument("-O", "--clobber", dest = "clobber", action = 'store_true', default = False, help = "Overwrite existing file if necessary.")

    parser.add_argument("--out-format", dest = "outformat", default = "NETCDF4_CLASSIC", help = "File output format (e.g., NETCDF3_CLASSIC, NETCDF4_CLASSIC, NETCDF4;pncgen only)", type = str, choices = 'NETCDF3_CLASSIC NETCDF4_CLASSIC NETCDF4'.split() + _writernames)

    parser.add_argument("--mode", dest = "mode", type = str, default = "w", help = "File mode for writing (w, a or r+ or with unbuffered writes ws, as, or r+s; pncgen only).", choices = 'w a r+ ws as r+s'.split())

    parser.add_argument('outpath', default = None, type = str, help='path to a output file formatted as --out-format')

def add_interactive_options(parser):
    try:
        parser.add_argument("-i", "--interactive", dest = "interactive", action = 'store_true', default = False, help = "Use interactive mode")
    except:
        pass

def add_plot_options(parser):
    add_interactive_options(parser)
    parser.add_argument("--matplotlibrc", dest = "matplotlibrc", type = lambda x: x.split(';'), action = 'append', default = [], help = 'rc options for matplotlib')
    parser.add_argument("--figure-keywords", dest = 'figure_keywords', type = lambda x: eval('dict(' + x + ')'), default = dict(), help = 'options for figure')
    
    parser.add_argument("--axes-keywords", dest = 'axes_keywords', type = lambda x: eval('dict(' + x + ')'), default = dict(), help = 'options for axes')
    
    parser.add_argument("--plot-commands", dest = "plotcommands", type = str, action = 'append', default = [], help = "Plotting functions to call for all variables expressions to execute in the context of the file. The figure object will always be 'fig'.")

    parser.add_argument("--figformat", dest = "figformat", type = str, default = 'png', help = "Any format supported by matplotlib")

    parser.add_argument("--norm", dest = "normalize", type = str, default = None, help = "Typical examples Normalize(), LogNorm(), BoundaryNorm([0, 10, 20, 30, 40], ncolors = 256)")

    parser.add_argument("--colorbar-formatter", dest = "colorbarformatter", type = str, default = None, help = "Typical examples LogFormatter(labelOnlyBase = False), ScalarFormatter(), '%%3g'")
    
    parser.add_argument("--overlay", dest = "overlay", action = 'store_true', help = "Enable overlay by setting equal to True")

    parser.add_argument('--no-squeeze', dest = 'squeeze', default = True, action = 'store_false', help = 'Squeeze automatically removes singleton dimensions; disabling requires user to remove singleton dimensions with --remove-singleton option')
    
    parser.add_argument("--figroot", help = 'Path for base of figure (suffixed by action)')

_pncepilog = """
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

def getparser2(actions = False):
    parent_parser = PNCArgumentParser(description = "PseudoNetCDF", 
                            formatter_class = RawDescriptionHelpFormatter, add_help = False, prog = 'PNC')
    if actions:
        subs = add_action_commands(parent_parser, withbasic = True, add_help = False)
    else:
        add_basic_options(parent_parser)
    return parent_parser


def getparser(has_ofile, plot_options = False, map_options = False, interactive = False, actions = False):
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
    parser = PNCArgumentParser(description = """PseudoNetCDF Argument Parsing


""", formatter_class=RawDescriptionHelpFormatter)
    parser.epilog = _pncepilog
    add_basic_options(parser)
    if interactive:
        add_interactive_options(parser)
        
    if actions:
        add_action_commands(parser, add_help = False)
    else:
        # Only has output file if pncgen is called.
        if has_ofile:
            add_output_options(parser)
        else:
            add_dump_options(parser)

        if plot_options:
            add_plot_options(parser)
        if map_options:
            add_map_options(parser)

    return parser
        
        
def pncparse(has_ofile = False, plot_options = False, map_options = False, interactive = False, args = None, parser = None, ifiles = []):
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
    inputargs = args
    if parser is None:
        parser = getparser(has_ofile, plot_options = plot_options, map_options = map_options, interactive = interactive)
    subparser = getparser(has_ofile = False, plot_options = plot_options, interactive = interactive)
    try:
        args = parser.parse_args(args = args)
    except Exception as e:
        raise e
    except BaseException as be:
        return [], args
    
    args.inputargs = inputargs or sys.argv
    inpaths = getattr(args, 'ifiles', [])
    inpnc = getattr(args, 'pnc', [])
    inopts = len(inpaths + inpnc)
    if len(ifiles) == 0 and inopts == 0:
        if not getattr(args, 'help', False):
            parser.print_help()
        
        print('\n*ifile or pnc are required')
        return [], args
    
    subargs = split_positionals(subparser, args)
    ifiles = [ifile for ifile in ifiles]
    ipaths = ['unknown'] * len(ifiles)
    for subarg in subargs:
        ipaths.extend(subarg.ifiles)
        ifiles.extend(pncprep(subarg)[0])
    
    # ifile is set to results from parsing
    # this includes the standard ifile
    args.ifiles = ifiles
    args.ipath = ipaths
    #if args.stack is not None:
    #    pass
    #elif len(args.ifiles) == 0 or (len(args.ifiles) - len(args.operators) - has_ofile) > 1:
    #    print('Too many arguments', len(args.ifiles), len(args.operators), has_ofile)
    #    parser.print_help()
    #    return
    if has_ofile or not getattr(args, 'outpath', None) is None:
        if not args.clobber and os.path.exists(args.outpath):
            parser.error(message = 'Output path (%s) exists; enable clobber (--clobber or -O) to overwrite.' % (args.outpath,))

    if plot_options or hasattr(args, 'matplotlibrc'):
        from matplotlib import rcParams
        for rcassigns in args.matplotlibrc:
            for rcassign in rcassigns:
                key, value = rcassign.split('=')
                rcParams[key] = value

    out = pncprep(args)
    return out

def add_action_commands(parser, withbasic = False, add_help = True):
    subparsers = parser.add_subparsers(help='Adding commands for subparsing', dest = 'subcommand')
    dumpparser = subparsers.add_parser('dump', description = 'Output Common Data Language (like from ncdump)', help = 'sub-command help', add_help = add_help)
    add_dump_options(dumpparser)
    if withbasic:
        add_basic_options(dumpparser)
    
    genparser = subparsers.add_parser('gen', description = 'Like ncgen with more format support', help = 'sub-command help', add_help = add_help)
    if withbasic:
        add_basic_options(genparser)
    add_output_options(genparser)

    evalparser = subparsers.add_parser('eval', description = 'Evaluation options.', help = 'sub-command help', add_help = add_help)
    add_eval_options(evalparser)
    if withbasic:
        add_basic_options(evalparser)
    
    mapparser = subparsers.add_parser('map', description = 'Make maps quickly', help = 'sub-command help', add_help = add_help)
    add_plot_options(mapparser)
    add_map_options(mapparser)
    if withbasic:
        add_basic_options(mapparser)
    
    for plotcmd in _plotcmds:
        plotparser = subparsers.add_parser(plotcmd, description = 'Plot data', add_help = add_help)
        add_plot_options(plotparser)
        if withbasic:
            add_basic_options(plotparser)
        if plotcmd == 'plotprofile':
            from .plotutil import add_vertprofile_options
            add_vertprofile_options(plotparser)
    


def do_actions(outargs):
    if not 'subcommand' in outargs:
        return
    if outargs.subcommand == 'dump':
        from .pncdump import pncdump
        for ifile in outargs.ifiles:
            pncdump(ifile, header = outargs.header, full_indices = outargs.full_indices, line_length = outargs.line_length, float_precision = outargs.float_precision, double_precision = outargs.double_precision, timestring = outargs.timestring, name = outargs.cdlname)
    
    elif outargs.subcommand == 'gen':
        from .pncgen import pncgen
        if len(outargs.ifiles) != 1:
            raise IOError('pncgen can output only 1 file; user requested %d' % len(outargs.ifiles))
        ifile, = outargs.ifiles
        pncgen(ifile, outargs.outpath, outmode = outargs.mode, format = outargs.outformat, verbose = outargs.verbose)
    elif outargs.subcommand == 'eval':
        from .pnceval import pnceval
        if len(outargs.ifiles) != 2:
            raise IOError('pnceval requires 2 files; user requested %d' % len(outargs.ifiles))
        pnceval(outargs)
    elif outargs.subcommand == 'map':
        if getattr(outargs, 'outpath', None) is None:
            outargs.outpath = outargs.figroot
        from .plotutil.pncmap import makemaps
        makemaps(outargs)        
    elif outargs.subcommand in _plotcmds:
        if getattr(outargs, 'outpath', None) is None:
            outargs.outpath = outargs.figroot
        from . import plotutil
        getattr(plotutil, outargs.subcommand)(outargs)        

def PNC(*args, **kwds):
    """
Arguments:
    args - Command Line arguments/options for PseudoNetCDF
           for full list of potential args PNC('--help')
    ifiles - (optional) pre-loaded input files
    actions -  (default: False)
        False: only open files do not make outputs
        True: enable dump,gen,map,etc output actions
              for action options see subparsers help
              e.g., PNC('dump', '--help', actions = True)

Returns:
    out - Namespace object with parsed arguments
          including a list of processed files (out.ifiles)

Example:
    # Single File
    out = PNC('--format=netcdf', inpath)    
    infile = out.ifiles[0]
    O3 = infile.variables['O3']

    # Multiple Files
    out = PNC('--format=netcdf', inpath1, inpath2)    
    infile1, infile2 = out.ifiles
    O3_1 = infile1.variables['O3']
    O3_2 = infile2.variables['O3']

    # With Actions
    out = PNC('dump', '--variables=O3', '--format=netcdf', inpath)    
PseudoNetCDF.icarttfiles.ffi1001.ffi1001 icartt/dc3-mrg60-dc8_merge_20120518_R7_thru20120622.ict {
dimensions:
        POINTS = 6817 ;

variables:
        double O3(TSTEP, LAY, ROW, COL);
                O3_ESRL:units = "ppbV" ;
                O3_ESRL:standard_name = "Ozone" ;
                O3_ESRL:missing_value = -999999 ;
    ...
    """
    ifiles = kwds.pop('ifiles', [])
    actions = kwds.pop('actions', None) 
    if actions is None:
        actions = args[0] in _allcmds

    parser = getparser2(actions = actions)
    nfiles = len(ifiles)
    ifiles, outargs = pncparse(args = args, ifiles = ifiles, parser = parser, **kwds)
    if nfiles > 0 and len(ifiles) != nfiles:
        raise IOError('PNC can only use files or paths, not both at this time')
    do_actions(outargs)
    
    return outargs

def pnc(*args, **kwds):
    """
Arguments - see PNC
Returns:
    file(s) - single file or, if more than 1 file is returned a list of files
    """
    ifiles = kwds.pop('ifiles', [])
    actions = kwds.pop('actions', None) 
    out = PNC(*args, ifiles = ifiles, actions = actions, **kwds).ifiles
    if len(out) == 1:
       return out[0]
    else:
       return out

def split_positionals(parser, args):
    import shlex
    positionals = args.ifiles
    parser.set_defaults(**dict([(k, v) for k, v in args._get_kwargs() if args.inherit or k == 'format']))
    ins = [shlex.split(pnc) for pnc in args.pnc]
    last_split = 0
    for i in range(len(positionals)):
        if positionals[i] in ('--', '--sep'):
            ins.append(positionals[last_split:i])
            last_split = i + 1
    ins.append(positionals[last_split:])
    try:
        outs = list(map(parser.parse_args, ins))
        for outns, inargs in zip(outs, ins):
            outns.inputargs = inargs
    except Exception as e:
        print(str(e))
        raise e
    except BaseException as be:
        print(str(be))
        return [], args
    
    for out in outs:
        if out.cdlname is None:
            out.cdlname = ', '.join(out.ifiles)
    if len(outs) == 1 and getattr(args, 'cdlname', None):
        args.cdlname = outs[0].cdlname
    return outs


def pncprep(args):
    #nifiles = len(args.ifiles) - has_ofile
    ipaths = args.ifiles[:]
    fs = getfiles(ipaths, args)
    fs = subsetfiles(fs, args)
    fs = seqpncbo(args.operators, fs, coordkeys = args.coordkeys)
    for expr in args.expressions:
        fs = [pncexpr(expr, f) for f in fs]
    for script in args.expressionscripts:
        expr = open(script).read()
        fs = [pncexpr(expr, f) for f in fs]
    args.ifiles = fs
    if getattr(args, 'cdlname', None) is None:
        try:
            args.cdlname = args.ipath[0]
        except:
            args.cdlname = 'unknown'
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
        if len(args.extractfile) > 0:
            extractfiles = []
            import shlex
            for extractfile in args.extractfile:
                extractfile, options = pncparse(has_ofile = False, args = shlex.split(extractfile))
                extractfiles.extend(extractfile)
            
            f = extract_from_file(f, extractfiles, method = args.extractmethod)
        for rd in args.removesingleton:
            f = removesingleton(f, rd)
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
        elif isinstance(ipath, (str,)) :
            try:
                allreaders = getreaderdict()
                if file_format in allreaders:
                    f = allreaders[file_format](ipath, **format_options)
                else:
                    f = eval(file_format)(ipath, **format_options)
            except Exception as e:
                oute = IOError('Unable to open path with %s(path, **%s)\n\tpath="%s"\n\terror="%s"' % (file_format, str(format_options), ipath, str(e)))
                raise oute# from e
        else:
            warn('File is type %s, which is unknown' % type(ipath))
            f = ipath
        
        history = getattr(f, 'history', getattr(f, 'HISTORY', ''))
        history += ' '.join(args.inputargs) + ';'
        laddconv = args.fromconv is not None and args.toconv is not None
        lslice = len(args.slice + args.reduce) > 0
        lexpr = len(args.expressions) > 0
        if args.variables is not None:
            f = getvarpnc(f, args.variables, coordkeys = args.coordkeys)
        #elif laddconv or lslice or lexpr:
        #    f = getvarpnc(f, None)
        for opts in args.attribute:
            add_attr(f, opts)
        for opts in args.masks:
            f = mask_vals(f, opts, metakeys = args.coordkeys)
        if laddconv:
            f = WrapPNC(f)
            try:
                eval('add_%s_from_%s' % (args.toconv, args.fromconv))(f, coordkeys = args.coordkeys)
            except Exception as e:
                warn('Cannot add %s from %s; %s' % (args.toconv, args.fromconv, str(e)))
        
        try:
            setattr(f, 'history', history)
        except: pass
        if args.mangle:
            f = manglenames(f)
        for rename in args.rename:
            f = pncrename(f, rename)
        fs.append(f)
    if args.stack is not None:
        fs = [stack_files(fs, args.stack, coordkeys = args.coordkeys)]
    if args.merge:
        fs = [merge(fs)]
    return fs
