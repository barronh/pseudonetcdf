__all__=['humidity']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- humidity Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              humidity files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from PseudoNetCDF.camxfiles.one3d.Memmap import one3d
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]
    
class humidity(one3d):
    """
    humidity provides a PseudoNetCDF interface for CAMx
    humidity files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> humidity_path = 'camx_humidity.bin'
        >>> rows,cols = 65,83
        >>> humidityfile = humidity(humidity_path,rows,cols)
        >>> humidityfile.variables.keys()
        ['TFLAG', 'HUM']
        >>> v = humidityfile.variables['HUM']
        >>> tflag = humidityfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> humidityfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    var_name='HUM'
    units='ppm'

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testHUM(self):
        import PseudoNetCDF.testcase
        humfile=humidity(PseudoNetCDF.testcase.CAMxHumidity,65,83)
        humfile.variables['TFLAG']
        self.assert_((humfile.variables['HUM'].mean(0).mean(1).mean(1)==array([ 28327.7109375 ,  28029.41601562,  27792.1484375 ,  27532.6171875 ,
         27164.37109375,  26810.52929688,  26387.734375  ,  25897.83984375,
         25292.76171875,  24546.08789062,  23261.65039062,  21335.16210938,
         19408.87109375,  17686.9296875 ,  15984.23535156,  14368.6796875 ,
         12869.12011719,  11707.51171875,  10261.96582031,   9015.51757812,
         8091.34423828,   7053.94873047,   5459.20556641,   3705.63745117,
         1855.64489746,    525.72961426,     94.93149567,     29.73270798],dtype='f')).all())

               
if __name__ == '__main__':
    unittest.main()
