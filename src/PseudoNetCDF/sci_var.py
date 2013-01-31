__doc__ = r"""
Scientific data has dimensions that have physical meaning and values
only have meaning in the context of their units.  This module implements
numpy arrays that are aware of their dimensions trying to vaguely adhere
to the Common Data Model from Unitdata at UCAR.

Each variable has as a property its dimensions names (dimensions).  Further,
each dimension name exists as a property and contains a one dimensional array
of values associated with that dimension.

For the purposes of ease of use, the standard properties of netCDF files
are attached and the arrays implement the Scientific.IO.NetCDF.NetCDFVariable
interfaces.
"""

__all__ = ['PseudoNetCDFFile', 'PseudoNetCDFVariableConvertUnit', 'PseudoNetCDFFileMemmap', 'PseudoNetCDFVariable', 'PseudoNetCDFMaskedVariable', 'PseudoIOAPIVariable', 'PseudoNetCDFVariables', 'Pseudo2NetCDF']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from numpy import array, zeros, ndarray, isscalar
from numpy.ma import MaskedArray, zeros as mazeros
from collections import defaultdict
from warnings import warn
    
from numpy import arange
from tempfile import NamedTemporaryFile as tnf
from types import MethodType
from units import convert
import operator,re,tempfile,warnings,sys,unittest

class PseudoNetCDFDimension(object):
    """
    Dimension object responds like that of netcdf4-python
    """
    def __init__(self, group, name, size):
        self._len = int(size)
    def isunlimitied(self):
        return False
    def __len__(self):
        return self._len
        
def PseudoNetCDFVariableConvertUnit(var,outunit):
    """
    Convert the unit of var and update the 
    associated IOAPI metadata
    """
    do = PseudoNetCDFFile()
    shape=var.shape
    for i,d in enumerate(var.dimensions):
        do.createDimension(d, shape[i])
    outvar=PseudoNetCDFVariable(do,var.long_name.strip(),var.typecode(),var.dimensions,values=convert(var,var.units,outunit))
    for k,v in var.__dict__.iteritems():
        setattr(outvar,k,v)
    outvar.units=outunit
    return outvar
    
class PseudoNetCDFFile(object):
    """
    PseudoNetCDFFile provides an interface and standard set of
    methods that a file should present to act like a netCDF file
    using the Scientific.IO.NetCDF.NetCDFFile interface.
    """
    def __new__(cls, *args, **kwds):
        new = object.__new__(cls, *args, **kwds)
        new.variables={}
        new.dimensions={}
        new._ncattrs = ()
        return new
    
    def __init__(self, *args, **properties):
        for k, v in properties.iteritems():
            setattr(self, k, v)

    def __setattr__(self, k, v):
        if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
            self._ncattrs += (k,)
        object.__setattr__(self, k, v)
    def createDimension(self,name,length):
        """
        name - string name for dimension
        length - maximum length of dimension
        """
        self.dimensions[name]=PseudoNetCDFDimension(self, name, length)

    def createVariable(self, name, type, dimensions, **properties):
        """
        name - string
        type - numpy dtype code (e.g., 'f', 'i', 'd')
        dimensions - tuple of dimension keys that can be
                     found in objects' dimensions dictionary
        """
        var = self.variables[name] = PseudoNetCDFVariable(self,name,type,dimensions, **properties)
        return var

    def close(self):
        """
        Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
        """
        pass

    def ncattrs(self):
        return self._ncattrs

    sync=close
    flush=close

class PseudoNetCDFFileMemmap(PseudoNetCDFFile):
    """
    Provides basic PseudoNetCDFFile functionality, but
    does not require that variables be created in memmory
    """
    def createVariable(self,name,type,dimensions,map,keep=True):
        var=PseudoNetCDFVariableMemmap(self,name,type,dimensions,map)
        if keep:
            self.variables[name]=var
        return var

class PseudoNetCDFVariable(ndarray):
    """
    PseudoNetCDFVariable presents the Scientific.IO.NetCDF.NetCDFVariable interface,
    but unlike that type, provides a contructor for variables that could be used
    without adding it to the parent file
    """
    def __setattr__(self, k, v):
        if not hasattr(self, k) and k[:1] != '_':
            self._ncattrs += (k,)
        ndarray.__setattr__(self, k, v)
    def ncattrs(self):
        return self._ncattrs
    def __new__(subtype,parent,name,typecode,dimensions,**kwds):
        """
        Creates a variable using the dimensions as defined in
        the parent object

        parent: an object with a dimensions variable
        name: name for variable
        typecode: numpy style typecode
        dimensions: a typle of dimension names to be used from
                    parrent
        kwds: Dictionary of keywords to be added as properties
              to the variable.  **The keyword 'values' is a special
              case that will be used as the starting values of
              the array

        """
        if 'values' in kwds.keys():
            result=kwds.pop('values')
        else:
            shape=[]
            for d in dimensions:
                dim = parent.dimensions[d]

                # Adding support for netCDF3 dimension objects
                if not isinstance(dim, int):
                    dim = len(dim)
                shape.append(dim)

            result=zeros(shape,typecode)

        result=result[...].view(subtype)

        if hasattr(result, '__dict__'):
            result.__dict__['typecode'] = lambda: typecode
            result.__dict__['dimensions'] = tuple(dimensions)
        else:
            result.__dict__={
                'typecode': lambda: typecode,
                'dimensions': tuple(dimensions)
            }

#        object.__setattr__(result, '_ncattrs', ())

        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        assert(hasattr(self, '_ncattrs') == False)
        self._ncattrs = ()
        if obj is None: return
        if hasattr(obj, '_ncattrs'):
            for k in obj._ncattrs:
                setattr(self, k, getattr(obj, k))
        
    
    def getValue(self):
        """
        Return scalar value
        """
        return self.item()

    def assignValue(self,value):
        """
        assign value to scalar variable
        """
        self.itemset(value)

class PseudoNetCDFMaskedVariable(MaskedArray, PseudoNetCDFVariable):
    def __new__(subtype,parent,name,typecode,dimensions,**kwds):
        """
        Creates a variable using the dimensions as defined in
        the parent object

        parent: an object with a dimensions variable
        name: name for variable
        typecode: numpy style typecode
        dimensions: a typle of dimension names to be used from
                    parrent
        kwds: Dictionary of keywords to be added as properties
              to the variable.  **The keyword 'values' is a special
              case that will be used as the starting values of
              the array

        """
        if 'values' in kwds.keys():
            result=kwds.pop('values')
        else:
            shape=[]
            for d in dimensions:
                dim = parent.dimensions[d]

                # Adding support for netCDF3 dimension objects
                if not isinstance(dim, int):
                    dim = len(dim)
                shape.append(dim)

            result=mazeros(shape,typecode)

        result=result[...].view(subtype)

        if hasattr(result, '__dict__'):
            result.__dict__['typecode'] = lambda: typecode
            result.__dict__['dimensions'] = tuple(dimensions)
        else:
            result.__dict__={
                'typecode': lambda: typecode,
                'dimensions': tuple(dimensions)
            }

#        object.__setattr__(result, '_ncattrs', ())

        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        MaskedArray.__array_finalize__(self, obj)
        assert(hasattr(self, '_ncattrs') == False)
        self._ncattrs = ()
        if obj is None: return
        if hasattr(obj, '_ncattrs'):
            for k in obj._ncattrs:
                setattr(self, k, getattr(obj, k))
        

    
def PseudoIOAPIVariable(parent,name,typecode,dimensions,**kwds):
    """
    Creates a variable using the dimensions as defined in
    the parent object
    
    parent: an object with a dimensions variable
    name: name for variable
    typecode: numpy style typecode
    dimensions: a typle of dimension names to be used from
                parent
    units: default = none
    long_name: default = name
    var_desc: default = name
    """

    retval = PseudoNetCDFVariable(parent, name, typecode, dimensions, **kwds)

    if not kwds.has_key('units'):
        warn('IOAPI variables must have units; %s has been initialized with "None" units')
        retval.units = 'None'
        
    if not kwds.has_key('long_name'):
        retval.long_name = name.ljust(16)

    if not kwds.has_key('var_desc'):
        retval.var_desc = name.ljust(80)

    return retval
    
class PseudoNetCDFVariables(defaultdict):
    """
    PseudoNetCDFVariables provides a special implementation
    of the default dictionary that provides efficient access 
    to variables of a PseudoNetCDFFile.  PseudoNetCDFFiles may
    have large variables that should only be loaded if accessed.
    PseudoNetCDFVariables allows a user to specify a function
    that can create variables on demand.
    """
    def __init__(self,func,keys):
        """
        func: Function that takes a key and provides a 
              PseudoNetCDFVariable
        keys: list of keys that the dictionary should
              act as if it has
        """
        self.__func=func
        self.__keys=keys
    def __missing__(self,k):
        """
        If the dictionary does not have a key, check if the
        user has provided that key.  If so, call the user 
        specifie function to create the variable.
        """
        if k in self.keys():
            return self.__func(k)
        else:
            raise KeyError, 'missing "%s"' % (k,)

    def addkey(self,k):
        """
        Allow the user to extend keys after the object
        has been created.
        """
        self.__keys.append(k)

    def keys(self):
        return self.__keys
    
    def has_key(self,k):
        return k in self.keys()
    
    def iteritems(self):
        for k in self.keys():
            yield k,self[k]
    
    def iterkeys(self):
        for k in self.keys():
            yield k

def get_ncf_object(path_or_object, mode, format = 'NETCDF4'):
    from os.path import exists, isfile
    from netcdf import NetCDFFile
    read_only = ('r', 'r+', 'rs', 'rs+', 'r+s')
    if isinstance(path_or_object, str):
        if exists(path_or_object):
            if isfile(path_or_object):
                ncf_object = NetCDFFile(path_or_object, mode, format = format)
            elif isdir(path_or_object):
                raise ValueError, "Got directory at %s; not sure what to do" % path_or_object
            else:
                raise ValueError, "Expected file or directory at %s" % path_or_object
        elif mode not in read_only:
            ncf_object = NetCDFFile(path_or_object, mode, format = format)
        else:
            raise IOError, "Cannot open missing file for reading"
    elif isinstance(path_or_object, NetCDFFile) or isinstance(path_or_object, PseudoNetCDFFile):
        return path_or_object
    elif path_or_object is None and mode not in read_only:
        tfile=tnf(mode='w+b')
        npath=tfile.name
        ncf_object=NetCDFFile(npath,mode)
    else:
        raise ValueError, "Not a path; not a netCDF file; not a PseudoNetCDF file... I don't know what to do"
    return ncf_object

def get_dimension_length(pfile, key):
    dim = pfile.dimsensions[key]
    if dim is None:
        for k in pfile.variables.keys():
            v = pfile.variables[k]
            if key in v.dimensions:
                return v[:].shape[list(v.dimensions).index(key)]
        return 0
    elif isinstance(dim, int):
        return dim
    else:
        return len(d)

class Pseudo2NetCDF:
    """
    Pseudo2NetCDF is a base class for conversion.  Properties and methods can
    be overwritten to facilitate conversion of special PseudoNetCDFFiles.
    
    Specifically: ignore_global_properties and ignore_variable_properties lists
    can be overwritten so that class properties and methods are not written
    to a netCDF file
    """
    ignore_global_re=re.compile('^_\w*(__\w*)?')
    ignore_variable_re=re.compile('^_\w*(__\w*)?')
    ignore_global_properties=['variables','dimensions']
    ignore_variable_properties=['typecode','dimensions']
    unlimited_dimensions = []
    create_variable_kwds = {}
    def convert(self,pfile,npath=None, inmode = 'r', outmode = 'w', format = 'NETCDF4'):
        pfile = get_ncf_object(pfile, inmode)
        nfile = get_ncf_object(npath, outmode, format = format)
        self.addDimensions(pfile,nfile)
        self.addGlobalProperties(pfile,nfile)
        self.addVariables(pfile,nfile)
        nfile.sync()
        return nfile
        
    def addDimensions(self,pfile,nfile):
        for d,v in pfile.dimensions.iteritems():
            if d in self.unlimited_dimensions:
                v = None
            elif not isinstance(v, (int, long)) and v is not None:
                v = len(v)
            
            nfile.createDimension(d,v)
        nfile.sync()
    
    def addGlobalProperties(self,pfile,nfile):
        for k in [k for k in pfile.__dict__.keys() if k not in self.ignore_global_properties and self.ignore_global_re.match(k)==None]:
            value=getattr(pfile,k)
            if not isinstance(value, MethodType):
                try:
                    setattr(nfile,k,value)
                except:
                    warn("Could not add %s to file" % k)

    def addVariableProperties(self,pvar,nvar):
        for a in [k for k in pvar.__dict__.keys() if k not in self.ignore_variable_properties and self.ignore_variable_re.match(k)==None]:
            value=getattr(pvar,a)
            if not isinstance(value, MethodType):
                try:
                    setattr(nvar,a,value)
                except:
                    if 'long_name' in pvar.__dict__.keys():
                        warn("Could not add %s=%s to variable %s" % (a,str(value),str(pvar.long_name)))
                    else:
                        warn("Could not add %s=%s to variable" % (a,str(value)))
                        
    
    def addVariable(self,pfile,nfile,k):
        pvar=pfile.variables[k]
        try:
            typecode = pvar.typecode()
        except:
            typecode = pvar[:].dtype.char
            
        nvar=nfile.createVariable(k,typecode,pvar.dimensions, **self.create_variable_kwds)

        if isscalar(nvar):
            nvar.assignValue(pvar)
        elif isinstance(pvar[:], MaskedArray):
            nvar[:] = pvar[:].filled()
        else:
            nvar[:] = pvar[:]
            
        self.addVariableProperties(pvar,nvar)
        nfile.sync()
        try:
            nfile.flush()
        except:
            pass
        del pvar,nvar

    def addVariables(self,pfile,nfile):
        for k in pfile.variables.keys():
            self.addVariable(pfile,nfile,k)
        nfile.sync()

class PseudoNetCDFTest(unittest.TestCase):
    from numpy import arange
    def setUp(self):
        self.tncf=PseudoNetCDFFile()
        
    def testNetCDFFileInit(self):
        tncf=self.tncf
        self.assert_(tncf.variables=={})
        self.assert_(tncf.dimensions=={})

        tncf.createDimension('TIME',24)
        tncf.createDimension('LAY',4)
        tncf.createDimension('ROW',5)
        tncf.createDimension('COL',6)

        self.assert_(len(tncf.dimensions['TIME'])==24)
        self.assert_(len(tncf.dimensions['LAY'])==4)
        self.assert_(len(tncf.dimensions['ROW'])==5)
        self.assert_(len(tncf.dimensions['COL'])==6)
        
        tncf.fish=2
        setattr(tncf,'FROG-DOG','HAPPY')

        o3=tncf.createVariable('O3','f',('TIME','LAY','ROW','COL'))
        self.assert_(tncf.variables.keys()==['O3'])
        
        o3[:] = arange(24*4*5*6).reshape(24,4,5,6)
        o3.units='ppbv'
        self.assert_((o3==arange(24*4*5*6).reshape(24,4,5,6)).all())
        self.assert_((tncf.variables['O3']==arange(24*4*5*6).reshape(24,4,5,6)).all())
        
        self.assert_(o3.typecode()=='f')

        filedims=tncf.dimensions.keys()
        filedims.sort()
        vardims=list(o3.dimensions)
        vardims.sort()

        self.assert_(filedims==vardims)
        
        n=Pseudo2NetCDF().convert(tncf)
        self.assertEqual(n.variables.keys(),['O3'])
        self.assertEqual(dict([(k, len(v)) for k, v in n.dimensions.iteritems()]),{'TIME': 24, 'LAY': 4, 'ROW': 5, 'COL': 6})
        self.assert_((n.variables['O3'][:]==tncf.variables['O3'][:]).all())
        self.assert_(n.variables['O3'].units=='ppbv')
        self.assertEqual(n.fish,2)
        self.assertEqual(getattr(n,'FROG-DOG'),'HAPPY')
        
    def testNetCDFVariables(self):
        tncf=PseudoNetCDFFile()
        tncf.createDimension('TIME',24)
        tncf.createDimension('LAY',4)
        tncf.createDimension('ROW',5)
        tncf.createDimension('COL',6)

        const=lambda *args,**kwds: PseudoNetCDFVariable(tncf,args[0],'f',('TIME','LAY','ROW','COL'),values=arange(24*4*5*6).reshape((24,4,5,6)))
        tncf.variables=PseudoNetCDFVariables(const,['NO','O3'])
        self.assert_((tncf.variables['O3']==arange(24*4*5*6).reshape(24,4,5,6)).all())
        
        
    def runTest(self):
        pass

if __name__ == '__main__':
    unittest.main()
