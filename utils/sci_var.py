HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from numpy import array, zeros,ndarray
try:
    from collections import defaultdict
except:
    from util import defaultdict
try:
    from Scientific.IO.NetCDF import NetCDFFile as ncf
except:
    from pynetcdf import NetCDFFile as ncf
    
from numpy import arange
from tempfile import NamedTemporaryFile as tnf
from types import MethodType
from pyPA.utils.units import convert
from pyPA.utils.util import AttrDict
import operator,re,tempfile,warnings,sys,unittest

__doc__="""
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

def PseudoNetCDFVariableConvertUnit(var,outunit):
    do=AttrDict({'dimensions':{}})
    shape=var.shape
    for i,d in enumerate(var.dimensions):
        do.dimensions[d]=shape[i]
    outvar=PseudoNetCDFVariable(do,var.long_name.strip(),var.typecode(),var.dimensions,convert(var,var.units,outunit))
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
    def __init__(self):
        self.variables={}
        self.dimensions={}
    
    def createDimension(self,name,length):
        self.dimensions[name]=length
    
    def createVariable(self,name,type,dimensions,keep=True):
        var=PseudoNetCDFVariable(self,name,type,dimensions)
        if keep:
            self.variables[name]=var
        return var

    def close(self):
        pass

    sync=close
    flush=close

class PseudoNetCDFFileMemmap(PseudoNetCDFFile):
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
    def __new__(subtype,parent,name,typecode,dimensions,values=None):
        """
        Creates a variable using the dimensions as defined in
        the parent object
        
        parent: an object with a dimensions variable
        name: name for variable
        typecode: numpy style typecode
        dimensions: a typle of dimension names to be used from
                    parrent
        """
        shape=[]
        for d in dimensions:
          shape.append(parent.dimensions[d])
        
        if values==None:
        	result=zeros(shape,typecode).view(subtype)
        else:
        	result=values.view(subtype)

        result.__dict__={
            'typecode': lambda: typecode,
            'dimensions': tuple(dimensions),
          }
        return result

    def getValue(self):
        """
        Return a array style variable
        """
        return self.view(ndarray)
  
    def assignValue(self,value):
        """
        assign value to variable value
        """
        self[...]=value

class PseudoNetCDFVariables(defaultdict):
    def __init__(self,func,keys):
            self.__func=func
            self.__keys=keys
    def __missing__(self,k):
    		if k in self.keys():
	            return self.__func(k)
	        else:
	            raise KeyError
    def keys(self):
            return self.__keys

class Pseudo2NetCDF:
    """
    Pseudo2NetCDF is a base class for conversion.  Properties and methods can
    be overwritten to facilitate conversion of special PseudoNetCDFFiles.
    
    Specifically: ignore_global_properties and ignore_variable_properties lists
    can be overwritten so that class properties and methods are not written
    to a netCDF file
    """
    ignore_global_re=re.compile('^_\w*__\w*')
    ignore_variable_re=re.compile('^_\w*__\w*')
    ignore_global_properties=['variables','dimensions']
    ignore_variable_properties=['typecode','dimensions']
    def convert(self,pfile,npath=None):
        if npath==None:
            tfile=tnf(mode='w+b')
            npath=tfile.name

        nfile=ncf(npath,'w')
        self.addDimensions(pfile,nfile)
        self.addGlobalProperties(pfile,nfile)
        self.addVariables(pfile,nfile)
        return nfile
        
    def addDimensions(self,pfile,nfile):
        for d,v in pfile.dimensions.iteritems():
            if d=='TSTEP' and type(nfile)!=PseudoNetCDFFile:
                nfile.createDimension(d,None)
            else:
                nfile.createDimension(d,v)
        nfile.sync()
    
    def addGlobalProperties(self,pfile,nfile):
        for k in [k for k in pfile.__dict__.keys() if k not in self.ignore_global_properties and self.ignore_global_re.match(k)==None]:
            value=getattr(pfile,k)
            if type(value)!=MethodType:
                setattr(nfile,k,value)

    def addVariableProperties(self,pvar,nvar):
        for a in [k for k in pvar.__dict__.keys() if k not in self.ignore_variable_properties and self.ignore_variable_re.match(k)==None]:
            value=getattr(pvar,a)
            if type(value)!=MethodType:
                setattr(nvar,a,value)
    
    def addVariable(self,pfile,nfile,k):
        pvar=pfile.variables[k]
        nvar=nfile.createVariable(k,pvar.typecode(),pvar.dimensions)
        nvar.assignValue(pvar)
        self.addVariableProperties(pvar,nvar)
        nfile.sync()
        nfile.flush()
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

        self.assert_(tncf.dimensions['TIME']==24)
        self.assert_(tncf.dimensions['LAY']==4)
        self.assert_(tncf.dimensions['ROW']==5)
        self.assert_(tncf.dimensions['COL']==6)
        
        tncf.fish=2
        setattr(tncf,'FROG-DOG','HAPPY')

        o3=tncf.createVariable('O3','f',('TIME','LAY','ROW','COL'))
        self.assert_(tncf.variables.keys()==['O3'])
        
        o3.assignValue(arange(24*4*5*6).reshape(24,4,5,6))
        o3.units='ppbv'
        self.assert_((o3.getValue()==arange(24*4*5*6).reshape(24,4,5,6)).all())
        self.assert_((tncf.variables['O3'].getValue()==arange(24*4*5*6).reshape(24,4,5,6)).all())
        
        self.assert_(o3.typecode()=='f')

        filedims=tncf.dimensions.keys()
        filedims.sort()
        vardims=list(o3.dimensions)
        vardims.sort()

        self.assert_(filedims==vardims)
        
        n=Pseudo2NetCDF().convert(tncf)
        self.assertEqual(n.variables.keys(),['O3'])
        self.assertEqual(n.dimensions,{'TIME': 24, 'LAY': 4, 'ROW': 5, 'COL': 6})
        self.assert_((array(n.variables['O3'].getValue())==tncf.variables['O3'].getValue()).all())
        self.assert_(n.variables['O3'].units=='ppbv')
        self.assertEqual(n.fish[...],2)
        self.assertEqual(getattr(n,'FROG-DOG'),'HAPPY')
        
    def testNetCDFVariables(self):
        tncf=PseudoNetCDFFile()
        tncf.createDimension('TIME',24)
        tncf.createDimension('LAY',4)
        tncf.createDimension('ROW',5)
        tncf.createDimension('COL',6)

        const=lambda *args,**kwds: arange(reduce(operator.mul,args[1:])).reshape((kwds['t'],kwds['l'],kwds['r'],kwds['c']))
        tncf.variables=PseudoNetCDFVariables(tncf,['NO','O3'],('TIME','LAY','ROW','COL'),const,(24,4,5,6),{'t':24,'l':4,'r':5,'c':6},'f')
        self.assert_((tncf.variables['O3'].getValue()==arange(24*4*5*6).reshape(24,4,5,6)).all())
        
        
    def runTest(self):
        pass

if __name__ == '__main__':
    unittest.main()
