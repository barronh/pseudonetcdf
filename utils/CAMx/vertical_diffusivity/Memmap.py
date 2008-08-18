HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['vertical_diffusivity']
#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd
from pyPA.utils.FortranFileUtil import OpenRecordFile,Int2Asc
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.CAMx.one3d.Memmap import one3d
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]
    
class vertical_diffusivity(one3d):
    var_name='KV'
    units='m**2/s'

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
               
    def testKV(self):
        vdfile=vertical_diffusivity('../../../../testdata/met/camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke',65,83)
        vdfile.variables['TFLAG']
        self.assert_((vdfile.variables['KV'].mean(0).mean(1).mean(1)==array([  13.65080357,   34.39198303,   68.02783966,   95.5898819 ,
          109.25765991,  112.92014313,  108.32209778,   97.25794983,
          84.1328125 ,   65.92033386,   46.97774506,   25.8343792 ,
          9.80327034,    2.89653206,    1.26993668,    1.12098336,
          1.13557184,    1.13372564,    1.19559622,    1.1675849 ,
          1.18877947,    1.18713808,    1.02371764,    1.02544105,
          1.21638143,    1.34624374,    1.03213251,    1.        ],dtype='f')).all())
    
if __name__ == '__main__':
    unittest.main()
