__all__=['landuse']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- landuse Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              landuse files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]
    
class landuse(PseudoNetCDFFile):
    """
    landuse provides a PseudoNetCDF interface for CAMx
    landuse files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> landuse_path = 'camx_landuse.bin'
        >>> rows,cols = 65,83
        >>> landusefile = landuse(landuse_path,rows,cols)
        >>> landusefile.variables.keys()
        ['TFLAG', 'FLAND', 'TOPO']
        >>> tflag = landusefile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = landusefile.variables['FLAND']
        >>> v.dimensions
        ('LANDUSE', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> landusefile.dimensions
        {'LANDUSE': 11, 'ROW': 65, 'COL': 83, 'VAR': 2}
    """
    
    def __init__(self,rf,rows=None,cols=None):
        self.rffile=OpenRecordFile(rf)
        
        self.dimensions=dict(ROW=rows,COL=cols,LANDUSE=11)
        self.variables={}
        self.__addvars()
        
    def __addvars(self):
        fland=PseudoNetCDFVariable(self,'FLAND','f',('LANDUSE','ROW','COL'),values=memmap(self.rffile.infile.name,'>f','r',0,(self.dimensions['ROW']*self.dimensions['COL']*self.dimensions['LANDUSE'],)).reshape((self.dimensions['LANDUSE'],self.dimensions['ROW'],self.dimensions['COL'])))
        fland.units='none'
        fland.long_name='FLAND'.ljust(16)
        fland.var_desc=fland.long_name
        self.variables['FLAND']=fland
        
        self.rffile.infile.seek(0,2)
        rflen=self.rffile.infile.tell()
        topo_start=fland.size*4+8
        if not rflen==topo_start:
            self.variables['TOPO']=PseudoNetCDFVariable(self,'TOPO','f',('ROW','COL'),values=memmap(self.rffile.infile.name,'>f','r',topo_start,(self.dimensions['ROW']*self.dimensions['COL'],)).reshape((self.dimensions['ROW'],self.dimensions['COL'])))
            topo=self.variables['TOPO']
            topo.units='none'
            topo.long_name='TOPO'.ljust(16)
            topo.var_desc=fland.long_name

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
               
    def testLU(self):
        import PseudoNetCDF.testcase
        lufile=landuse(PseudoNetCDF.testcase.CAMxLandUse,65,83)
        self.assert_((lufile.variables['FLAND'].mean(1).mean(1)==array([ 0.02158891223, 0.07436342537, 0.01858459227, 0.0509638451 , 0.08883236349, 0.17303463817, 0.33603012562, 0.00156044273, 0.05768220499, 0.1666303575 , 0.01072917785],dtype='f')).all())

if __name__ == '__main__':
    unittest.main()
