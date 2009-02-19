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

__all__=['pncdump','dump_from_cmd_line', 'pncdump_parser', 'pncdump_parseri']
from numpy import float32, float64, int32, int64, uint32, uint64, ndenumerate, savetxt
from warnings import warn
from optparse import OptionParser

import textwrap
import sys
import os

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
                   bool = "%s")
    float_fmt = "%%.%df" % (float_precision,)
    int_fmt = "%i"
    # initialize indentation as 8 characters
    # based on ncdump
    indent = 8*" "
    
    # Global and variable properties must be an instance 
    # of one of the following
    property_types = (str, int, float)
    
    # First line of CDL
    sys.stdout.write("netcdf %s {\n" % (name,))
    
    ###########################
    # CDL Section 1: dimensions
    ###########################
    sys.stdout.write("dimensions:\n")
    for dim_name, dim in f.dimensions.iteritems():
        sys.stdout.write(1*indent+("%s = %s ;\n" % (dim_name,dim)))
    
    ###################################
    # CDL Section 2: variables metadata
    ###################################
    sys.stdout.write("\nvariables:\n")
    for var_name, var in f.variables.iteritems():
        var_type = dict(float32='float', \
                        float64='double', \
                        int32='integer', \
                        int64='long', \
                        bool='bool')[var.dtype.name]
        sys.stdout.write(1*indent+("%s %s%s;\n" % (var_type, var_name,str(var.dimensions).replace('\'',''))))
        for prop_name, prop in var.__dict__.iteritems():
            if isinstance(prop,property_types):
                sys.stdout.write(2*indent+("%s:%s = \"%s\" ;\n" % (var_name,prop_name,prop)))
    
    ################################
    # CDL Section 3: global metadata
    ################################
    sys.stdout.write("\n\n// global properties:\n")
    for prop_name, prop in f.__dict__.iteritems():
        if isinstance(prop,property_types):
            sys.stdout.write(2*indent+(":%s = %s ;\n" % (prop_name, repr(prop))))

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
                        if i == (var.size-1):
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
                    for hri, hr in enumerate(var):
                        for layi, lay in enumerate(hr):
                            for rowi, row in enumerate(lay):
                                fmt = ', '.join(row.size * [formats[var.dtype.name]])
                                array_str = fmt % tuple(row)
                                if hri == (f.dimensions['TSTEP']-1) and \
                                   layi == (f.dimensions['LAY']-1) and \
                                   rowi == (f.dimensions['ROW']-1):
                                    array_str += ";"
                                else:
                                    array_str += ","

                                try:
                                    sys.stdout.write(textwrap.fill(array_str, line_length, initial_indent = '  ', subsequent_indent = '    '))
                                    sys.stdout.write('\n')
                                except IOError:
                                    exit()
        except KeyboardInterrupt:
            sys.stdout.flush()
            exit()
    sys.stdout.write("}\n")

def dump_from_cmd_line(file_path, options, reader = None):
    """
    dump_from_cmd_line provides basic error checking, but requires
    a reader function to initialize the PseudoNetCDFFile
    
    parser  - OptParse parser object with basic implementation 
              taken from pncdump_parser
    options - OptParse options object implementing, at minimum
              the pncdump_parser values
    args    - OptParse args list implementing, at minimum
              the pncdump_parser expected values
    reader  - Any PseudoNetCDFFile initializer that takes only
              args[0] as an argument and returns a PseudoNetCDFFile
    """
    pncdump(reader(file_path), name = options.name, header = options.header, variables = options.variables, line_length = options.line_length, full_indices = options.full_indices, float_precision = options.float_precision, double_precision = options.double_precision)

if __name__ == '__main__':
    from PseudoNetCDF.pncdump import pncdump_parser, \
                                    dump_from_cmd_line
    parser = pncdump_parser()
    parser.add_argument("cols", int)
    parser.add_argument("rows", int)
    from PseudoNetCDF.camxfiles.vertical_diffusivity.Memmap import vertical_diffusivity
    (file_path, options, extra_args_dict) = parser.parse_args()

    dump_from_cmd_line(file_path, options, lambda path: vertical_diffusivity(path, **extra_args_dict))
