__doc__ = r"""
.. _MetaNetCDF
:mod:`MetaNetCDF` -- PseudoNetCDF manipulation utilities
========================================================

.. module:: MetaNetCDF
   :platform: Unix, Windows
   :synopsis: Provides tools for manipulating NetCDF-like files
              and variables.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

__all__ = ['add_derived', 'time_avg_new_unit', 'window', 'newresolution', 'MetaNetCDF', 'WindowFromFile', 'file_master']

HeadURL = "$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from numpy import array, where, logical_or, repeat, mean, sum, zeros
#This Package modules
from PseudoNetCDF.sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariable, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariableConvertUnit
from PseudoNetCDF.ArrayTransforms import CenterTime
                            

class add_derived(PseudoNetCDFFile):
    """add_derived provides a simple interface to add derived variables
    to a PseudoNetCDFFile interface
    
    create a new class with the following modifications:
    overwrite __childclass__ with the base class
    overwrite __addvars__ with a list of keys for variable names you 
                          intend to derive
    for each key name, create a key interface funciton (e.g. key = DEPTH, interface = __DEPTH__)
    """
    __childclass__ = None
    __addvars__ = {}
    def __init__( * args, **kwds):
        self = args[0]
        self.__child = self.__childclass__( * args[1:], **kwds)
        self.dimensions = self.__child.dimensions
        self.variables = PseudoNetCDFVariables(self.__variables__, self.__child.variables.keys() + self.__addvars__)
    def __variables__(self, key):
        if key in self.__addvars__:
           return getattr(self, '__' + key + '__')()
        else:
            return self.__child.variables[key]

class time_avg_new_unit(PseudoNetCDFFile):
    """
    This base class provides an interface for converting
    instantaneous data to time averaged.  It also provides
    an optional unit conversion.
    
    ex:
        >>> windfile = wind(wind_path, 65, 83)
        >>> windfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
        >>> windfile.variables['U'].units
        'm/s'
        >>> windfile.variables['V'].units
        'm/s'
        >>> class windavgnewunit(time_avg_new_unit):
        ...     __reader__ = wind
        >>> windavgfile = windavgnewunit(wind_path, rows, cols, {'U': 'km/h', 'V': 'mph'})
        >>> windavgfile.dimensions
        {'TSTEP': 24, 'LAY': 28, 'ROW': 65, 'COL': 83}
        >>> windfile.variables['U'].units
        'km/h'
        >>> windfile.variables['V'].units
        'mph'
    """
    __reader__ = None
    def __init__(self, rffile, rows, cols, outunit = {}, endhour = True):
        self.__file = self.__reader__(rffile, rows, cols)
        self.createDimension('TSTEP', len(self.__file.dimensions['TSTEP'] - 1))
        self.createDimension('LAY', len(self.__file.dimensions['LAY']))
        self.createDimension('ROW', len(self.__file.dimensions['ROW']))
        self.createDimension('COL', len(self.__file.dimensions['COL']))
        self.createDimension('VAR', len(self.__file.dimensions['VAR']))
        if self.__file.dimensions.has('DATE-TIME'):
            self.createDimension('DATE-TIME', len(self.__file.dimensions['DATE-TIME']))
        else:
            self.createDimension('DATE-TIME', 2)
        self.__outunit = outunit
        self.variables = PseudoNetCDFVariables(self.__variables, self.__file.variables.keys())
        self.__timeslice = {True:slice(1, None), False:slice(None, -1)}[endhour]
        v = self.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'), keep = True)
        v[:] = self.__file.variables['TFLAG'][self.__timeslice]
        v.long_name = 'Time flag'
        v.units = 'DATE-TIME'

    def __variables(self, k):
        outunit = self.__outunit.get(k, None)
        var = self.__file.variables[k]
        if outunit==None:
            outunit = var.units
        return PseudoNetCDFVariableConvertUnit(self.__decorator(var, PseudoNetCDFVariable(self, k, var.typecode(), var.dimensions, values = CenterTime(var))), outunit)
    
    def __decorator(self, ovar, nvar):
        for a, v in ovar.__dict__.iteritems():
            setattr(nvar, a, v)
        return nvar

class window(PseudoNetCDFFile):
    """
    Window takes a netcdf or PseudoNetCDF file and creates
    a windowed version.
    
    ex:
        >>> windfile = wind(wind_path, 65, 83)
        >>> subsetfile = window(windfile, tslice = slice(8, 18), \
                                          kslice = slice(0, 1), \
                                          jslice = slice(19, 37), \
                                          islice = slice(19, 37))
        >>> windfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
        >>> windfile.variables['U'].shape
        (25, 28, 65, 83)
        >>> subsetfile.dimensions
        {'TSTEP': 10, 'LAY': 1, 'ROW': 18, 'COL': 18}
        >>> windfile.variables['U'].shape
        (10, 1, 18, 18)
        
    """
    def __init__(self, ncffile, tslice = slice(None), kslice = slice(None), jslice = slice(None), islice = slice(None)):
        self.__file = ncffile
        self.__idx = [tslice, kslice, jslice, islice]
        self.dimensions = self.__file.dimensions.copy()
        any_non_time_key = [k for k in self.__file.variables.keys() if 'TFLAG' not in k][0]
        self.dimensions['TSTEP'], self.dimensions['LAY'], self.dimensions['ROW'], self.dimensions['COL'] \
                            = map(lambda x: PseudoNetCDFDimension(None,  None, x), self.__file.variables[any_non_time_key][self.__idx].shape)
        self.variables = PseudoNetCDFVariables(self.__variables, self.__file.variables.keys())
        
    def __variables(self, k):
        if 'TFLAG' in k:
            return self.__file.variables[k]
        
        ov = self.__file.variables[k]
        nv = ov[self.__idx]
        from PseudoNetCDF.sci_var import Pseudo2NetCDF
        Pseudo2NetCDF().addVariableProperties(nv, ov)
        return nv

class newresolution(PseudoNetCDFFile):
    """
    newresolution converts a netcdf or PseudoNetCDF file
    to a new resolution based on user input.  This class
    updates dimensions and variables, but metadata is unaffected.
    
    ex:
        >>> kvfile = kv(kv_path, 65, 83)
        >>> newresfile = newresolution(kvfile, (2, 3), 4000, 1000)
        >>> kvfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
        >>> kvfile.variables['KV'].shape
        (25, 28, 65, 83)
        >>> newresfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 18, 'COL': 18}
        >>> newresfile.variables['KV'].shape
        (25, 28, 260, 332)
    """
    def __init__(self, ncffile, dimension, oldres, newres, repeat_method = repeat, condense_method = sum, nthick = 0):
        PseudoNetCDFFile.__init__(self)
        self.__dimension = array(dimension, ndmin = 1)
        oldres = array(oldres, ndmin = 1)
        newres = array(newres, ndmin = 1)
        self.__mesh = newres / oldres.astype('f')
        self.__condense = condense_method
        self.__repeat = repeat_method
        self.__file = ncffile
        self.__nthick = nthick
            
        if not logical_or((self.__mesh % 1) == 0, (1. / self.__mesh) % 1 ==0).any():
            raise ValueError, "One resolution must be a factor of the other."

        Pseudo2NetCDF().addDimensions(self.__file, self)
        any_non_time_key = [k for k in self.__file.variables.keys() if 'TFLAG' not in k][0]
        for dk, dfactor in zip(self.__dimension, 1./self.__mesh):
            dimo = self.dimensions[dk]
            ndimo = self.createDimension(str(dk), len(dimo)*dfactor)
            ndimo.setunlimited(dimo.isunlimited())
        v = self.__file.variables[any_non_time_key]
        v = self.__method(v)
        
        self.variables = PseudoNetCDFVariables(self.__variables, self.__file.variables.keys())
 
    def __method(self, a):
        dimensions = self.__dimension
        meshs = self.__mesh
        for dimension, mesh in zip(dimensions, meshs):
            axisi = list(a.dimensions).index(dimension)
            newshape = list(a.shape)
            if mesh < 1:
                method = lambda a: self.__repeat(a, (1./mesh), axisi)
                newshape[axisi] = newshape[axisi] * mesh
            else:
                newshape[axisi:axisi + 1] = int(newshape[axisi] // mesh), int(mesh)
                method = lambda a: self.__condense(a.reshape(newshape), axisi + 1)
            a = method(a)
            
        return a

    def __variables(self, k):
        if 'TFLAG' in k and (self.__axis != 0).any():
            raise KeyError, "Tflag is off limits"
        else:
            ov = self.__file.variables[k]
            v = self.__method(ov)
            Pseudo2NetCDF().addVariableProperties(ov, v)
            return v
                
class MetaNetCDF(PseudoNetCDFFile):
    """
    MetaNetCDF provides a basic interface for combining files
    and adding derived variables.
    
    ex:
        >>> kvfile = kv(kv_path, 65, 83)
        >>> zpfile = height_pressure(zp_path, 65, 83)
        >>> windfile = wind(wind_path, 65, 83)
        >>> metfile = MetaNetCDF([kvfile, zpfile, windfile])
        >>> wind_speed_calc = lambda self: (self.variables['U'][:]**2 + \
                                            self.variables['V'][:]**2)**.5
        >>> metfile.addMetaVariable('WS', wind_speed_calc)
        >>> metfile.variables['WS'].shape
        (25, 28, 65, 83)
    """
    __metavars__ = {}
    def addMetaVariable(self, key, func):
        self.variables.addkey(key)
        self.__metavars__[key] = func

    def __init__(self, files):
        self.__files = files
        keys = []
        for f in self.__files:
            for k, d in f.dimensions.iteritems():
                if len(d)==1 and k=='LAY':
                    k = 'SURFLAY'
                if k not in self.dimensions.keys():
                    self.createDimension(k, len(d))
            keys.extend(f.variables.keys())
            for k in f.__dict__.keys():
                if k not in self.__dict__.keys():
                    setattr(self, k, getattr(f, k))
        keys = list(set(keys))
        self.variables = PseudoNetCDFVariables(self.__variables, keys)
        self.createDimension('VAR', len(keys) - 1)
    
    def __getattribute__(self, k):
        try:
            return PseudoNetCDFFile.__getattribute__(self, k)
        except AttributeError:
            for f in self.__files:
                try:
                    return getattr(f, k)
                except:
                    pass
            raise AttributeError, "%s not found" % k

    def childvariables(self, k):
        for f in self.__files:
            if k in f.variables.keys():
                v = f.variables[k]
                if k=='TFLAG':
                    v = PseudoNetCDFVariable(self, 'TFLAG', 'i', v.dimensions, values = v[:][:, [0], :].repeat(len(self.dimensions['VAR']), 1))
                    v.long_name = 'TFLAG'.ljust(16)
                    v.var_desc = 'TFLAG'.ljust(16)
                    v.units = 'DATE-TIME'
                    
                if k=='LAY' and k in self.dimensions.keys() and len(k.shape) > 1:
                    if v.shape[1]==1:
                        dims = list(v.dimensions)
                        dims[1] = 'SURFLAY'
                        v.dimensions = tuple(dims)
                return v
    
    def __variables(self, k):
        if k in self.__metavars__.keys():
            return self.__metavars__[k](self)
        else:
            return self.childvariables(k)
        raise KeyError, '%s not in any files' % k
file_master = MetaNetCDF

def WindowFromFile(WindowThis, WindowFrom):
    """
    WindowFromFile creates a bounding box to window one file
    from the meta data of another file.
    
    ex:
        >>> wind12k = wind(wind12k_path, 89, 89)
        >>> wind04k = wind(wind04k_path, 65, 83)
        >>> window4k_from_wind12k = WindowFromFile(wind12k, wind04k)
        >>> V12k = wind12k.variables['V']
        >>> V04k = wind04k.variables['V']
        >>> V04k_12k = window4k_from_wind12k.variables['V']
        >>> V04k_12k.shape
        (25, 28, 89, 89)
        >>> V04k_12k.shape
        (25, 28, 65, 83)
        >>> V04k_12k.shape # Note the trimmed buffer cells
        (25, 28, 21, 27)
    """
    xoffset = abs(WindowThis.XORIG - WindowFrom.XORIG)
    ioffset = 0
    while xoffset != 0 and  not xoffset % WindowThis.XCELL:
        xoffset += WindowFrom.XCELL
        ioffset += 1
    istart= xoffset / WindowThis.XCELL

    yoffset = abs(WindowThis.YORIG - WindowFrom.YORIG)
    joffset = 0
    while yoffset != 0 and not yoffset % WindowThis.YCELL:
        yoffset += WindowFrom.YCELL
        joffset += 1
    jstart= yoffset / WindowThis.YCELL

    iend = istart + int((WindowFrom.NCOLS - 2 * ioffset) * WindowFrom.XCELL / WindowThis.XCELL)
    jend = jstart + int((WindowFrom.NROWS - 2 * joffset) * WindowFrom.YCELL / WindowThis.YCELL)
    islice = slice(istart, iend)
    jslice = slice(jstart, jend)

    outf = MetaNetCDF([WindowThis])
    def AddMetaVar(k, jslice, islice):
        outf.addMetaVariable(k, lambda self: self.childvariables(k)[:, :, jslice, islice])
    for k in WindowThis.variables.keys():
        AddMetaVar(k, jslice, islice)
    return outf