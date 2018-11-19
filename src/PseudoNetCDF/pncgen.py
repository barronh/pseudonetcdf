from __future__ import print_function
import sys
from warnings import warn
from types import MethodType
from PseudoNetCDF.netcdf import NetCDFVariable
from .sci_var import PseudoNetCDFFile
from .sci_var import get_ncf_object
import numpy as np
from PseudoNetCDF.camxfiles import Writers as CAMxWriters
import PseudoNetCDF.geoschemfiles as geoschemwriters
import PseudoNetCDF.icarttfiles.ffi1001 as icarttwriters


if sys.version_info > (3,):
    long = int


class Pseudo2NetCDF:
    """
    Pseudo2NetCDF is a base class for conversion.  Properties and methods can
    be overwritten to facilitate conversion of special PseudoNetCDFFiles.

    Specifically: ignore_global_properties and ignore_variable_properties lists
    can be overwritten so that class properties and methods are not written
    to a netCDF file
    """
    import re
    ignore_global_re = re.compile(r'^_\w*(__\w*)?')
    ignore_variable_re = re.compile(r'^_\w*(__\w*)?')
    ignore_global_properties = ['variables', 'dimensions']
    ignore_variable_properties = ['typecode', 'dimensions']
    special_properties = ['_fillvalue', '_FillValue']
    unlimited_dimensions = []
    create_variable_kwds = {}

    def __init__(self, datafirst=False, verbose=1):
        self.datafirst = datafirst
        self.verbose = verbose

    def convert(self, pfile, npath=None, inmode='r', outmode='w',
                format='NETCDF4'):
        pfile = get_ncf_object(pfile, inmode)
        nfile = get_ncf_object(npath, outmode, format=format)
        if self.verbose:
            print("Adding dimensions", file=sys.stdout)
        self.addDimensions(pfile, nfile)
        if self.verbose:
            print("Adding globals", file=sys.stdout)
        self.addGlobalProperties(pfile, nfile)
        if self.verbose:
            print("Adding variables", file=sys.stdout)
        self.addVariables(pfile, nfile)
        nfile.sync()
        return nfile

    def addDimensions(self, pfile, nfile):
        for d in pfile.dimensions.keys():
            self.addDimension(pfile, nfile, d)

    def addDimension(self, pfile, nfile, d):
        v = pfile.dimensions[d]
        unlim = (d in self.unlimited_dimensions or v.isunlimited())
        if not isinstance(v, (int, long)) and v is not None:
            v = len(v)

        if unlim:
            if isinstance(nfile, PseudoNetCDFFile):
                nd = nfile.createDimension(d, v)
                nd.setunlimited(True)
            else:
                nd = nfile.createDimension(d, None)
        else:
            nd = nfile.createDimension(d, v)

        nfile.sync()

    def addGlobalProperties(self, pfile, nfile):
        for k in [k for k in pfile.ncattrs()
                  if (k not in self.ignore_global_properties and
                      self.ignore_global_re.match(k) is None)]:
            value = getattr(pfile, k)
            if not isinstance(value, MethodType):
                try:
                    setattr(nfile, k, value)
                except TypeError as e:
                    if isinstance(value, bool):
                        setattr(nfile, k, np.int8(value))
                    else:
                        raise e
                except Exception as e:
                    warn("Could not add %s to file; %s: %s" % (k, type(e), e))

    def addVariableProperties(self, pvar, nvar):
        for a in [k for k in pvar.ncattrs()
                  if ((k not in self.ignore_variable_properties and
                       self.ignore_variable_re.match(k) is None) or
                      k in self.special_properties)]:
            value = getattr(pvar, a)
            if isinstance(nvar, NetCDFVariable) and a == '_FillValue':
                continue
            if not isinstance(value, MethodType):
                try:
                    nvar.setncattr(a, value)
                    # setattr(nvar,a,value)
                except TypeError as e:
                    if isinstance(value, bool):
                        nvar.setncattr(a, np.int8(value))
                    else:
                        raise e
                except Exception as e:
                    if 'long_name' in pvar.ncattrs():
                        warn("Could not add %s=%s to variable %s; %s" %
                             (a, str(value), str(pvar.long_name), e))
                    else:
                        warn("Could not add %s=%s to variable; %s" %
                             (a, str(value), e))

    def addVariable(self, pfile, nfile, k, data=True):
        pvar = pfile.variables[k]
        try:
            typecode = pvar.typecode()
        except Exception:
            typecode = pvar[...].dtype.char

        create_variable_kwds = self.create_variable_kwds.copy()
        if hasattr(pvar, 'missing_value'):
            create_variable_kwds['fill_value'] = pvar.missing_value
        elif hasattr(pvar, 'fill_value'):
            create_variable_kwds['fill_value'] = pvar.fill_value
        elif hasattr(pvar, '_FillValue'):
            create_variable_kwds['fill_value'] = pvar._FillValue

        nvar = nfile.createVariable(
            k, typecode, pvar.dimensions, **create_variable_kwds)
        self.addVariableProperties(pvar, nvar)
        if data:
            self.addVariableData(pfile, nfile, k)
            nfile.sync()

        try:
            nfile.flush()
        except Exception:
            pass
        del pvar, nvar

    def addVariableData(self, pfile, nfile, k):
        from numpy.ma import MaskedArray
        from numpy import isscalar
        nvar = nfile.variables[k]
        pvar = pfile.variables[k]
        if isscalar(nvar) or nvar.ndim == 0:
            if isinstance(pvar, NetCDFVariable):
                pvar = pvar[...]
            nvar[...] = pvar
        elif isinstance(pvar[...], MaskedArray):
            nvar[:] = pvar[...].filled(getattr(nvar, 'fill_value', getattr(
                nvar, '_FillValue', getattr(pvar, 'missing_value', -9999))))
        else:
            nvar[:] = pvar[...]

    def addVariables(self, pfile, nfile):
        for k in pfile.variables.keys():
            if self.verbose:
                print("Defining", k, file=sys.stdout)
            self.addVariable(pfile, nfile, k, data=self.datafirst)

        if self.datafirst:
            nfile.sync()

        if not self.datafirst:
            for k in pfile.variables.keys():
                if self.verbose:
                    print("Populating", k, file=sys.stdout)
                self.addVariableData(pfile, nfile, k)
            nfile.sync()


def pywriter(ifile, outpath, data=True):
    print("""# Import Libraries and Functions
from netCDF4 import Dataset
from numpy import *
from numpy.ma import masked_array

# Open file for writing
outpath = '%s'
outfile = Dataset(outpath, 'w')

""" % outpath)
    print("# Define Dimensions")
    for dk, dv in ifile.dimensions.items():
        dl = len(dv)
        if dv.isunlimited():
            dl = None
        print("dim_%s = outfile.createDimension('%s', %s); # %d" %
              (dk, dk, dl, len(dv)))

    print("# Add global properties")
    for pk in ifile.ncattrs():
        pv = getattr(ifile, pk)
        if isinstance(pv, bool):
            pv = int(pv)
        print("setattr(outfile, '%s', %s)" % (pk, repr(pv)))

    print("## Define Variables")
    print("vars = {}")
    for vk, v in ifile.variables.items():
        print("# Defining " + vk)
        print("var = vars['%s'] = outfile.createVariable('%s', '%s', %s)" % (
            vk, vk, v.dtype.char, v.dimensions))
        for pk in v.ncattrs():
            pv = getattr(v, pk)
            print("setattr(var, '%s', %s)" % (pk, repr(pv)))
        print("")

    if data:
        print("## Populate Variables")
        for vk, v in ifile.variables.items():
            print("# Populating " + vk)
            if isinstance(v, np.ma.MaskedArray):
                vtype = np.ma.MaskedArray
            else:
                vtype = np.ndarray
            print("var = vars['%s']" % vk)
            print("var[:] = %s" % (repr(v[:].view(type=vtype))))


def pncgen(ifile, outpath, inmode='r', outmode='w', format='NETCDF4_CLASSIC',
           verbose=1, complevel=0, writer_kw=None):
    """
    ifile - input file to write out
    outpath - path to outputfile
    inmode - how is file read (if ifile is a path)
    outmode - w, w+s
    format - any PseudoNetCDF or Dataset option
    """
    if writer_kw is None:
        writer_kw = {}
    if format[:6] == 'NETCDF':
        p2n = Pseudo2NetCDF()
        p2n.verbose = verbose
        if complevel > 0:
            p2n.create_variable_kwds['zlib'] = True
            p2n.create_variable_kwds['complevel'] = complevel
        return p2n.convert(ifile, outpath, inmode=inmode, outmode=outmode,
                           format=format)

    from ._getwriter import getwriterdict
    writerdict = getwriterdict()
    if format == 'python':
        pywriter(ifile, outpath)
    elif format == 'csv':
        from .textfiles._delimited import ncf2csv
        ncf2csv(ifile, outpath)
    elif format in writerdict:
        writer = writerdict[format]
        return writer(ifile, outpath, **writer_kw)
    else:
        for writers in [CAMxWriters, geoschemwriters, icarttwriters]:
            writer = getattr(writers, 'ncf2%s' % format, None)
            if writer is not None:
                return writer(ifile, outpath, **writer_kw)
                break
        else:
            raise KeyError('Unknown output file type "%s"' % format)


def main():
    from .pncparse import pncparse
    ifiles, options = pncparse(has_ofile=True, interactive=False)
    if len(ifiles) != 1:
        raise IOError(
            'pncgen can output only 1 file; user requested %d' % len(ifiles))
    ifile, = ifiles
    return pncgen(ifile, options.outpath, outmode=options.mode,
                  format=options.outformat, verbose=options.verbose), options


if __name__ == '__main__':
    main()
