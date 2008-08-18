HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxRead.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
from types import GeneratorType
import unittest
import struct,sys,os,operator
from warnings import warn
from tempfile import TemporaryFile as tempfile
import os,sys

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype
from pynetcdf import NetCDFFile as ncf

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd,timerange
from pyPA.utils.util import cartesian,sliceit
from pyPA.utils.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.CAMx.one3d.Read import one3d as one3d

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

__all__=['humidity']

class humidity(one3d):
    var_name='HUM'
    units='ppm'


class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testHUM(self):
        humfile=humidity('../../../../testdata/met/camx_hum.20000825.hgbpa_04km.TCEQuh1_eta.v43',65,83)
        self.assert_((humfile.variables['HUM'].mean(0).mean(1).mean(1)==array([ 28327.7109375 ,  28029.41601562,  27792.1484375 ,  27532.6171875 ,
         27164.37109375,  26810.52929688,  26387.734375  ,  25897.83984375,
         25292.76171875,  24546.08789062,  23261.65039062,  21335.16210938,
         19408.87109375,  17686.9296875 ,  15984.23535156,  14368.6796875 ,
         12869.12011719,  11707.51171875,  10261.96582031,   9015.51757812,
         8091.34423828,   7053.94873047,   5459.20556641,   3705.63745117,
         1855.64489746,    525.72961426,     94.93149567,     29.73270798],dtype='f')).all())
        
if __name__ == '__main__':
    unittest.main()
