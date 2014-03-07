import sys
from warnings import warn
from types import MethodType
from PseudoNetCDF.netcdf import NetCDFFile
from sci_var import PseudoNetCDFFile
from sci_var import get_ncf_object
import numpy as np

class Pseudo2NetCDF:
    """
    Pseudo2NetCDF is a base class for conversion.  Properties and methods can
    be overwritten to facilitate conversion of special PseudoNetCDFFiles.
    
    Specifically: ignore_global_properties and ignore_variable_properties lists
    can be overwritten so that class properties and methods are not written
    to a netCDF file
    """
    import re
    ignore_global_re=re.compile('^_\w*(__\w*)?')
    ignore_variable_re=re.compile('^_\w*(__\w*)?')
    ignore_global_properties=['variables','dimensions']
    ignore_variable_properties=['typecode','dimensions']
    special_properties = ['_fillvalue', '_FillValue']
    unlimited_dimensions = []
    create_variable_kwds = {}
    def __init__(self, datafirst = False, verbose = True):
        self.datafirst = datafirst
        self.verbose = verbose
    def convert(self,pfile,npath=None, inmode = 'r', outmode = 'w', format = 'NETCDF4'):
        pfile = get_ncf_object(pfile, inmode)
        nfile = get_ncf_object(npath, outmode, format = format)
        if self.verbose: print >> sys.stdout, "Adding dimensions"
        self.addDimensions(pfile,nfile)
        if self.verbose: print >> sys.stdout, "Adding globals"
        self.addGlobalProperties(pfile,nfile)
        if self.verbose: print >> sys.stdout, "Adding variables"
        self.addVariables(pfile,nfile)
        nfile.sync()
        return nfile
        
    def addDimensions(self,pfile,nfile):
        for d,v in pfile.dimensions.iteritems():
            unlim = (d in self.unlimited_dimensions or v.isunlimited())
            if not isinstance(v, (int, long)) and v is not None:
                v = len(v)
            
            if unlim:
                if isinstance(nfile, PseudoNetCDFFile):
                    nd = nfile.createDimension(d,v)
                    nd.setunlimited(True)
                else:
                    nd = nfile.createDimension(d,None)
            else:
                nd = nfile.createDimension(d,v)

        nfile.sync()
    
    def addGlobalProperties(self,pfile,nfile):
        for k in [k for k in pfile.ncattrs() if k not in self.ignore_global_properties and self.ignore_global_re.match(k)==None]:
            value=getattr(pfile,k)
            if not isinstance(value, MethodType):
                try:
                    setattr(nfile,k,value)
                except TypeError, e:
                    if isinstance(value, bool):
                        setattr(nfile, k, np.int8(value))
                    else:
                        raise e
                except Exception, e:
                    warn("Could not add %s to file; %s: %s" % (k, type(e), e))

    def addVariableProperties(self,pvar,nvar):
        for a in [k for k in pvar.ncattrs() if (k not in self.ignore_variable_properties and self.ignore_variable_re.match(k)==None) or k in self.special_properties]:
            value=getattr(pvar,a)
            if not isinstance(value, MethodType):
                try:
                    setattr(nvar,a,value)
                except TypeError, e:
                    if isinstance(value, bool):
                        setattr(nvar, a, np.int8(value))
                    else:
                        raise e
                except Exception, e:
                    if 'long_name' in pvar.ncattrs():
                        warn("Could not add %s=%s to variable %s; %s" % (a,str(value),str(pvar.long_name), e))
                    else:
                        warn("Could not add %s=%s to variable; %s" % (a,str(value), e))
                        
    
    def addVariable(self,pfile,nfile,k, data = True):
        pvar=pfile.variables[k]
        try:
            typecode = pvar.typecode()
        except:
            typecode = pvar[...].dtype.char
        
        create_variable_kwds = self.create_variable_kwds.copy()
        if hasattr(pvar, 'fill_value'):
            create_variable_kwds['fill_value'] = pvar.fill_value
        nvar=nfile.createVariable(k,typecode,pvar.dimensions, **create_variable_kwds)
        self.addVariableProperties(pvar,nvar)
        if data:
            self.addVariableData(pfile, nfile, k)
        nfile.sync()
        try:
            nfile.flush()
        except:
            pass
        del pvar,nvar

    def addVariableData(self, pfile, nfile, k):
        from numpy.ma import MaskedArray
        from numpy import ndarray, isscalar
        nvar = nfile.variables[k]
        pvar = pfile.variables[k]
        if isscalar(nvar) or nvar.ndim == 0:
            nvar.assignValue(pvar)
        elif isinstance(pvar[...], MaskedArray):
            nvar[:] = pvar[...].filled()
        else:
            nvar[:] = pvar[...]
        

    def addVariables(self,pfile,nfile):
        for k in pfile.variables.keys():
            if self.verbose: print >> sys.stdout, "Defining", k
            self.addVariable(pfile,nfile,k, data = self.datafirst)
        nfile.sync()
        if not self.datafirst:
            for k in pfile.variables.keys():
                if self.verbose: print >> sys.stdout, "Populating", k
                self.addVariableData(pfile,nfile,k)
            nfile.sync()

from PseudoNetCDF.camxfiles import Writers as CAMxWriters
import  PseudoNetCDF.geoschemfiles as geoschemwriters

def pncgen(ifile,outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4_CLASSIC', verbose = True):
    if format[:6] == 'NETCDF':
        p2n = Pseudo2NetCDF()
        p2n.verbose = verbose
        return p2n.convert(ifile, outpath, inmode = inmode, outmode = outmode, format = format)
    else:
        for writers in [CAMxWriters, geoschemwriters]:
            writer = getattr(writers, 'ncf2%s' % format, None)
            if not writer is None:
                return writer(ifile, outpath)
                break
        else:
            raise KeyError('Unknown output file type "%s"' % format)
            
            
        
    
def main():
    from pncparse import pncparser
    ifile, ofile, options = pncparser(has_ofile = True)
    return pncgen(ifile, ofile, outmode = options.mode, format = options.outformat, verbose = options.verbose), options

if __name__ == '__main__':
    main()
