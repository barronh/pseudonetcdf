__doc__ = r"""
.. _dumper
:mod:`dumper` -- PseudoNetCDF dump module
============================================

.. module:: dumper
   :platform: Unix, Windows
   :synopsis: Provides ncdump like functionaility for PseudoNetCDF
.. moduleauthor:: Barron Henderson <barronh@unc.edu>

Example implementation:
    if __name__ == '__main__':
        from PseudoNetCDF.pncdump import pncdump_parser, \
                                        dump_from_cmd_line
        parser = pncdump_parser()
        parser.add_argument("cols", int)
        parser.add_argument("rows", int)
        from PseudoNetCDF.camxfiles.vertical_diffusivity.Memmap import vertical_diffusivity
        (file_path, options, extra_args_dict) = parser.parse_args()
    
        dump_from_cmd_line(file_path, options, lambda path: vertical_diffusivity(path, **extra_args_dict))

Example implementation usage:
    python -m PseudoNetCDF.pncdump camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke 65 83
"""

__all__=['pncdump','dump_from_cmd_line', 'pncdump_parser']
from numpy import float32, float64, int16, int32, int64, uint32, uint64, ndenumerate, savetxt
from warnings import warn
from optparse import OptionParser

import textwrap
import sys
import os
import operator

class pncdump_parser(OptionParser):
    def __init__(self, *args, **kwds):
        OptionParser.__init__(self, *args, **kwds)
        self.set_conflict_handler("resolve")
        self.set_usage("Usage: pncdump [-h] [-v var1[,...]] [-f c|f]] [-l len] [-n name] [--help] file_path")
        self.add_option("-h", "--header", dest="header",action = "store_true", default=False)
        self.add_option("-v", "--variables", dest="variables",default="", metavar = "var[,...]")
        self.add_option("-n", "--name", dest="name",default=None, metavar = "filename")
        self.add_option("-l", "--length", dest="line_length",default=80, metavar = "len")
        self.add_option("-f", "--full", dest="full_indices",default=None, metavar = "[C|F]")
        self.add_option("-p", "--precision", dest="precision",default=None, metavar = "fdig,pdig")
        self.__extra_arguments = []
        
    def check_values(self, options, args):
        args_total = 1+len(self.__extra_arguments)
        if len(args) != args_total:
            self.error(msg="Requires file path pluss any arguments necessary initialize the file: %d" % (args_total,))
        else:
            if os.path.exists(args[0]):
                file_path = args[0]
            else:
                self.error(msg="The file_path argument must be a real file")

        extra_args = dict([(name,type_func(args[i+1])) for i,(name,type_func) in enumerate(self.__extra_arguments)])
            
    
        if options.variables is not None:
            options.variables = options.variables.split()
        
        if options.name is None:
            options.name = os.path.basename(args[0])
            
        options.line_length = int(options.line_length)
        
        if options.full_indices is not None:
            options.full_indices = (options.full_indices+'c')[0].lower()
        
        if options.precision is None:
            options.float_precision = 8
            options.double_precision = 16
        else:
            options.float_precision, options.double_precision = [int(v) for v in options.precision.split(',')]
        
        return file_path, options, extra_args
    
    def add_argument(self, name, type_func):
        self.__extra_arguments.append((name, type_func))
        self.set_usage(self.usage+ " " + name)
        

def pncdump(f, name = 'unknown', header = False, variables = [], line_length = 80, full_indices = None, float_precision = 8, double_precision = 16):
    """
    pncdump is designed to implement basic functionality
    of the NetCDF ncdump binary.
    
    f         - a PseudoNetCDFFile object
    name      - string name for the file 
                (equivalent to ncdump -n name)
    header    - boolean value for display of header only
                (equivalent to ncdump -h)
    variables - iterable of variable names for subsetting
                data display (equivalent to ncddump -v var[,...]
    """
    formats = dict(float64 = "%%.%de" % (double_precision,), \
                   float32 = "%%.%de" % (float_precision,), \
                   int32 = "%i", \
                   int64 = "%i", \
                   str = "%s", \
                   bool = "%s", \
                   string8 = "'%s'")
    float_fmt = "%%.%df" % (float_precision,)
    int_fmt = "%i"
    # initialize indentation as 8 characters
    # based on ncdump
    indent = 8*" "
    
    # Global and variable properties must be an instance 
    # of one of the following
    property_types = (str, int, int16, int32, int64, float, float32, float64)
    
    # First line of CDL
    sys.stdout.write("netcdf %s {\n" % (name,))
    
    ###########################
    # CDL Section 1: dimensions
    ###########################
    sys.stdout.write("dimensions:\n")
    for dim_name, dim in f.dimensions.iteritems():
        sys.stdout.write(1*indent+("%s = %s ;\n" % (dim_name,len(dim))))
    
    ###################################
    # CDL Section 2: variables metadata
    ###################################
    sys.stdout.write("\nvariables:\n")
    for var_name, var in f.variables.iteritems():
        var_type = dict(float32='float', \
                        float64='double', \
                        int32='integer', \
                        int64='long', \
                        bool='bool', \
                        string8='char')[var.dtype.name]
        sys.stdout.write(1*indent+("%s %s%s;\n" % (var_type, var_name,str(var.dimensions).replace('\'','').replace(',)',')'))))
        for prop_name, prop in var.__dict__.iteritems():
            if isinstance(prop,property_types):
                sys.stdout.write(2*indent+("%s:%s = \"%s\" ;\n" % (var_name,prop_name,prop)))
    
    ################################
    # CDL Section 3: global metadata
    ################################
    sys.stdout.write("\n\n// global properties:\n")
    for prop_name, prop in f.__dict__.iteritems():
        if isinstance(prop,property_types):
            sys.stdout.write(2*indent+(":%s = %s ;\n" % (prop_name, repr(prop).replace("'",'"'))))

    if not header:
        # Error trapping prevents the user from getting an error
        # when they cancel a dump or when they break a redirected
        # pipe
        try:
            #####################
            # CDL Section 4: data
            #####################
            sys.stdout.write("\n\ndata:\n")
            
            # data indentation is only 1 space
            indent = " "
            
            # Subset variables for output
            display_variables = [var_name for var_name in f.variables.keys() if var_name in variables or variables == []]
            if variables != []:
                if len(variables) < len(display_variables):
                    warn("Not all specified variables were available")
            
            # For each variable outptu data 
            # currently assumes 3-D data
            for var_name in display_variables:
                var = f.variables[var_name]
                sys.stdout.write(1*indent+("%s =\n" % var_name))
                if full_indices is not None:
                    id_display = {'f': lambda idx: str(tuple([idx[i]+1 for i in range(len(idx)-1,-1,-1)])), \
                                  'c': lambda idx: str(idx)}[full_indices]
                                  
                    for i, val in ndenumerate(var):
                        fmt = 2*indent+formats[var.dtype.name]
                        array_str = fmt % val
                        if i == tuple(map(lambda x_: x_ - 1, var.shape)):
                            array_str += ";"
                        else:
                            array_str += ","

                        array_str += " // %s%s \n" % (var_name, id_display(i))
                        try:
                            sys.stdout.write(array_str)
                        except IOError:
                            sys.stdout.close()
                            exit()
                else:
                    dimensions = [len(f.dimensions[d]) for d in var.dimensions]
                    if len(dimensions) > 1:
                        first_dim = reduce(operator.mul,dimensions[:-1])
                        second_dim = dimensions[-1]
                        shape = [first_dim, second_dim]
                    else:
                        shape = [1]+dimensions
                    var2d = var[...].reshape(shape)
                    fmt = ', '.join(shape[-1] * [formats[var.dtype.name]])
                    for rowi, row in enumerate(var2d):
                        try:
                            row = tuple(row)
                        except:
                            pass
                        array_str = fmt % row
                        if rowi == (shape[0]-1):
                            array_str += ';'
                        else:
                            array_str += ','
                            
                        try:
                            sys.stdout.write(textwrap.fill(array_str, line_length, initial_indent = '  ', subsequent_indent = '    '))
                            sys.stdout.write('\n')
                        except IOError, e:
                            warn(repr(e) + "; Typically from CTRL+C or exiting less")
                            exit()
                                            
                    
        except KeyboardInterrupt:
            sys.stdout.flush()
            exit()
    sys.stdout.write("}\n")

if __name__ == '__main__':
    from optparse import OptionParser
    from camxfiles.Memmaps import *
    from camxfiles.Readers import irr as irr_read, ipr as ipr_read
    from icarttfiles.ffi1001 import ffi1001
    from geoschemfiles import *
    try:
        from netCDF4 import Dataset as netcdf
    except:
        pass
    from sci_var import reduce_dim, slice_dim, getvarpnc, extract

    parser = OptionParser()
    parser.set_usage("""Usage: python -m %prog [-f uamiv|bpch|ffi1001|...] ifile [ofile]

    ifile - path to a file formatted as type -f
    ofile - path to the desired output
    
    -f --format - format of the file either uamiv (CAMx), bpch (GEOS-Chem) or ffi1001 (optionally has comma delimited arguments for opening the format)
    """)

    parser.add_option("-f", "--format", dest = "format", default = 'uamiv', help = "File format")
    
    parser.add_option("-v", "--variables", dest = "variables", default = None,
                        help = "Variable names or regular expressions (using match) separated by ','. If a group(s) has been specified, only variables in that (those) group(s) will be selected.")

    parser.add_option("-s", "--slice", dest = "slice", type = "string", action = "append", default = [],
                        help = "bpch variables have dimensions (time, layer, lat, lon), which can be subset using dim,start,stop,stride (e.g., --slice=layer,0,47,5 would sample every fifth layer starting at 0)")

    parser.add_option("-e", "--extract", dest = "extract", action = "append", default = [],
                        help = "lon/lat coordinates to extract")
    
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
    elif len(options.slice + options.reduce) > 0:
        f = getvarpnc(f, None)
    for opts in options.slice:
        f = slice_dim(f, opts)
    for opts in options.reduce:
        f = reduce_dim(f, opts)
    if len(options.extract) > 0:
        extract(f, options.extract)

    pncdump(f)
