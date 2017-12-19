__doc__ = r"""
.. _dumper
:mod:`dumper` -- PseudoNetCDF dump module
============================================

.. module:: pncdump
   :platform: Unix, Windows
   :synopsis: Provides ncdump like functionaility for PseudoNetCDF
.. moduleauthor:: Barron Henderson <barronh@unc.edu>

"""

__all__=['pncdump',]
from numpy import float32, float64, int16, int32, int64, ndindex, ndenumerate, nan, savetxt, isscalar, ndarray, ma, prod, set_printoptions, get_printoptions, array2string, inf
from warnings import warn
from collections import defaultdict
from PseudoNetCDF import PseudoNetCDFMaskedVariable

import sys
if sys.version_info.major == 3:
    from io import BytesIO as StringIO
    commaspace = u', '
    semicolon = b';'
else:
    from StringIO import StringIO
    commaspace = ', '
    semicolon = ';'
    BrokenPipeError = None
import textwrap
import operator

def exception_handler(e, outfile):
    """
    Exceptions are largely ignored
    """
    try:
        outfile.flush()
        outfile.close()
    except:
        pass
    if isinstance(e, KeyboardInterrupt) or isinstance(e, BrokenPipeError) or (isinstance(e, IOError) and e.strerror == 'Broken pipe'):
        exit()
    else:
        warn(repr(e) + "; Typically from CTRL+C or exiting less")
        exit()

def pncdump(f, name = 'unknown', header = False, variables = [], line_length = 80, full_indices = None, float_precision = 8, double_precision = 16, isgroup = False, timestring = False, outfile = sys.stdout):
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

    pncdump(vertical_diffusivity('camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke',rows=65,cols=83))
    """
    file_type = str(type(f)).split("'")[1]
    float_fmt = "%%.%dg" % (float_precision,)
    double_fmt = "%%.%dg" % (double_precision,)
    int_fmt = "%i"
    formats = defaultdict(lambda: "%s",
                   float = double_fmt, \
                   float64 = double_fmt, \
                   float32 = float_fmt, \
                   int32 = "%i", \
                   uint32 = "%i", \
                   int64 = "%i", \
                   str = "%s", \
                   bool = "%s", \
                   string8 = "'%s'")

    funcs = dict()#float = lambda x: double_fmt % x)    # initialize indentation as 8 characters
    # based on ncdump
    indent = 8*" "
    if isgroup:
        startindent = 4*" "
    else:
        startindent = 4*""
        
    # First line of CDL
    if not isgroup: outfile.write("%s %s {\n" % (file_type, name,))
    
    ###########################
    # CDL Section 1: dimensions
    ###########################
    outfile.write(startindent + "dimensions:\n")
    for dim_name, dim in f.dimensions.items():
        if dim.isunlimited():
            outfile.write(startindent + 1*indent+("%s = UNLIMITED // (%s currently) \n" % (dim_name,len(dim))))
        else:
            outfile.write(startindent + 1*indent+("%s = %s ;\n" % (dim_name,len(dim))))
    
    ###################################
    # CDL Section 2: variables metadata
    ###################################
    if len(f.variables) > 0:
        outfile.write("\n" + startindent + "variables:\n")
    for var_name, var in f.variables.items():
        var_type = dict(float32='float', \
                        float64='double', \
                        int32='integer', \
                        uint32='integer', \
                        int64='long', \
                        bool='bool', \
                        string8='char', \
                        string80='char').get(var.dtype.name, var.dtype.name)
        outfile.write(startindent + 1*indent+("%s %s%s;\n" % (var_type, var_name,str(var.dimensions).replace('u\'', '').replace('\'','').replace(',)',')'))))
        for prop_name in var.ncattrs():
            prop = getattr(var, prop_name)
            outfile.write(startindent + 2*indent+("%s:%s = %s ;\n" % (var_name,prop_name,repr(prop).replace("'", '"'))))
    
    ################################
    # CDL Section 3: global metadata
    ################################
    outfile.write("\n\n// global properties:\n")
    for prop_name in f.ncattrs():
        prop = getattr(f, prop_name)
        outfile.write(startindent + 2*indent+(":%s = %s ;\n" % (prop_name, repr(prop).replace("'",'"'))))

    if hasattr(f, 'groups'):
        for group_name, group in f.groups.items():
            outfile.write(startindent + 'group %s:\n' % group_name)
            pncdump(group, name = name, header = header, variables = variables, line_length = line_length, full_indices = full_indices, float_precision = float_precision, double_precision = double_precision, isgroup = True)
    if not header:
        # Error trapping prevents the user from getting an error
        # when they cancel a dump or when they break a redirected
        # pipe
        try:
            #####################
            # CDL Section 4: data
            #####################
            outfile.write("\n\n" + startindent + "data:\n")
            
            # data indentation is only 1 space
            indent = " "
            
            # Subset variables for output
            display_variables = [var_name for var_name in f.variables.keys() if var_name in variables or variables == []]
            if variables != []:
                if len(variables) < len(display_variables):
                    warn("Not all specified variables were available")
            
            # For each variable output data 
            # currently assumes 3-D data
            for var_name in display_variables:
                var = f.variables[var_name]
                if isinstance(var, PseudoNetCDFMaskedVariable) or hasattr(var, '_FillValue'):
                    def writer(row, last):
                        if isscalar(row) or row.ndim == 0:
                            outfile.write(startindent + '  ' + str(row.filled().astype(ndarray)) + ';\n')
                            return
                        #old = get_printoptions()
                        #set_printoptions(threshold = inf, linewidth = line_length)
                        #tmpstr =  startindent + '    ' + array2string(row, separator = commaspace, formatter = funcs).replace('\n', '\n' + startindent + '    ')[1:-1].replace('--', '_')
                        #if last:
                        #    tmpstr += ';'
                        #else:
                        #    tmpstr += commaspace                        
                        #set_printoptions(**old)
                        #outfile.write(tmpstr)
                        #outfile.write('\n')
                        
                        tmpstr = StringIO()
                        if ma.getmaskarray(row).all():
                            tmpstr.write(b', '.join([b'_'] * row.size) + b', ')
                        else:
                            savetxt(tmpstr, ma.filled(row), fmt, delimiter = commaspace, newline = commaspace)
                        if last:
                            tmpstr.seek(-2, 1)
                            tmpstr.write(semicolon)
                        tmpstr.seek(0, 0)
                        tmpstr = tmpstr.read().decode('ASCII')
                        tmpstr = tmpstr.replace(fmt % getattr(row, 'fill_value', 0) + ',', '_,')
                        tmpstr = textwrap.fill(tmpstr, line_length, initial_indent = startindent + '  ', subsequent_indent = startindent + '    ')
                        try:
                            outfile.write(tmpstr)
                            outfile.write('\n')
                        except Exception as e:
                            exception_handler(e, outfile)
                else:
                    def writer(row, last):
                        if isscalar(row) or row.ndim == 0:
                            outfile.write(startindent + '  ' + str(row.astype(ndarray)) + ';\n')
                            return
                        old = get_printoptions()
                        set_printoptions(threshold = inf, linewidth = line_length)
                        tmpstr =  startindent + '    ' + array2string(row, separator = commaspace, formatter = funcs).replace('\n', '\n' + startindent + '    ')[1:-1]
                        #tmpstr = StringIO()
                        #savetxt(tmpstr, row, fmt, delimiter = commaspace, newline =commaspace)
                        if last:
                            tmpstr += ';'
                            #tmpstr.seek(-2, 1)
                            #tmpstr.write(semicolon)
                        else:
                            tmpstr += commaspace
                        #tmpstr.seek(0, 0)
                        #outfile.write(textwrap.fill(str(tmpstr.read()), line_length, initial_indent = startindent + '  ', subsequent_indent = startindent + '    '))
                        set_printoptions(**old)
                        try:
                            outfile.write(tmpstr)
                            outfile.write('\n')
                        except Exception as e:
                            exception_handler(e, outfile)
                        
                        
                outfile.write(startindent + 1*indent+("%s =\n" % var_name))
                if var_name in ('time', 'time_bounds') and timestring:
                    from PseudoNetCDF.coordutil import gettimes, gettimebnds
                    if var_name == 'time':
                        times = gettimes(f)
                    elif var_name == 'time_bounds':
                        times = gettimebnds(f)
                    
                    for i in ndindex(var.shape):
                        val = var[i]
                        if ma.is_masked(val):
                            array_str = '_'
                        else:
                            array_str = startindent + 2*indent + str(val)

                        if i == tuple(map(lambda x_: x_ - 1, var.shape)):
                            array_str += ";"
                        else:
                            array_str += ","

                        array_str += " // %s%s %s \n" % (var_name, i, str(times[i]))
                        try:
                            outfile.write(array_str)
                        except Exception as e:
                            exception_handler(e, outfile)
                elif full_indices is not None:
                    id_display = {'f': lambda idx: str(tuple([idx[i]+1 for i in range(len(idx)-1,-1,-1)])), \
                                  'c': lambda idx: str(idx)}[full_indices]
                                  
                    #for i, val in ndenumerate(var):
                    for i in ndindex(var.shape):
                        val = var[i]
                        if ma.is_masked(val):
                            array_str = '_'
                        else:
                            fmt = startindent + 2*indent+formats[var.dtype.name]
                            array_str = fmt % val

                        if i == tuple(map(lambda x_: x_ - 1, var.shape)):
                            array_str += ";"
                        else:
                            array_str += ","

                        array_str += " // %s%s \n" % (var_name, id_display(i))
                        try:
                            outfile.write(array_str)
                        except Exception as e:
                            exception_handler(e, outfile)
                else:
                    dimensions = [len(f.dimensions[d]) for d in var.dimensions]
                    if len(dimensions) > 1:
                        first_dim = prod(dimensions[:-1])
                        second_dim = dimensions[-1]
                        shape = [first_dim, second_dim]
                    else:
                        shape = [1]+dimensions
                    var2d = var[...].reshape(*shape)
                    fmt = ', '.join(shape[-1] * [formats[var.dtype.name]])
                    fmt = formats[var.dtype.name]
                    lastrow = var2d.shape[0] - 1
                    for rowi, row in enumerate(var2d):
                        try:
                            writer(row, rowi == lastrow)
                        except Exception as e:
                            exception_handler(e, outfile)
                                            
                    
        except Exception as e:
            exception_handler(e, outfile)
            
    outfile.write("}\n")
    return outfile

def main():
    from .pncparse import pncparse
    ifiles, options = pncparse(has_ofile = False)
    for ifile in ifiles:
        pncdump(ifile, header = options.header, full_indices = options.full_indices, line_length = options.line_length, float_precision = options.float_precision, double_precision = options.double_precision, timestring = options.timestring, name = options.cdlname)

if __name__ == '__main__':
    main()
