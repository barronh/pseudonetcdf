HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['cloud_rain']
#Distribution packages
import unittest
import struct
from warnings import warn

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd
from pyPA.utils.FortranFileUtil import OpenRecordFile,Int2Asc
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]
    
class cloud_rain(PseudoNetCDFFile):
    def __init__(self,rf,rows=None,cols=None):
        f=file(rf,'rb')
        f.seek(0,2)
        flen=f.tell()
        offset=struct.unpack('>i',file(rf,'r').read(4))[0]+8
        self.__memmap=memmap(rf,'>f','r',offset=offset)
        self.dimensions={}
        ncols,nrows,nlays=self.dimensions['COL'],self.dimensions['ROW'],self.dimensions['LAY']=struct.unpack({35:'>i15ciiii',40:'>i20ciiii'}[offset],file(rf,'r').read(offset))[-4:-1]
        self.STIME,self.SDATE=struct.unpack({35:'>i15ciiiiifi',40:'>i20ciiiiifi'}[offset],file(rf,'r').read(offset+12))[-2:]
        if self.SDATE<10000:
            self.SDATE+=2000000
        if (ncols!=cols and cols!=None) or (rows!=rows and rows!=None):
            warn('Files says cols=%d, rows=%d, and lays=%d; you said cols=%d and rows=%d' % (ncols,nrows,nlays,cols,rows))
            
        self.dimensions['DATE-TIME']=2
        self.VERSION,varkeys={35:('<4.3',['CLOUD','PRECIP','COD','TFLAG']),40:('4.3',['CLOUD','RAIN','SNOW','GRAUPEL','COD','TFLAG'])}[offset]
        self.dimensions['TSTEP']=int((flen-offset)/((len(varkeys)-1)*nlays*(nrows*ncols+2)*4+16))
        self.dimensions['VAR']=len(varkeys)-1
        
        self.NVARS=self.dimensions['VAR']
        self.NLAYS=self.dimensions['LAY']
        self.NROWS=self.dimensions['ROW']
        self.NCOLS=self.dimensions['COL']
        self.FTYPE=1
        
        self.variables=PseudoNetCDFVariables(self.__var_get,varkeys)
        
        self.SDATE,self.STIME=self.variables['TFLAG'][0,0,:]

    def __set_var(self,key,vals_idx):
        times=self.dimensions['TSTEP']
        lays=self.dimensions['LAY']
        rows=self.dimensions['ROW']
        cols=self.dimensions['COL']
        v=PseudoNetCDFVariable(self,key,'f',('TSTEP','LAY','ROW','COL'),self.__memmap[vals_idx].reshape(times,lays,rows,cols))
        v.units={'COD':'None'}.get(key,'g/m**3')
        v.long_name=key
        v.var_desc=key
        self.variables[key]=v
        
    def __var_get(self,key):
        times=self.dimensions['TSTEP']
        rows=self.dimensions['ROW']
        cols=self.dimensions['COL']
        lays=self.dimensions['LAY']
        vars=len(self.variables.keys())-1
        hour=1
        date=2
        cloud=3
        rain=4
        snow=5,
        graupel=6,
        cod=7
        stagger=8
        out_idx=zeros(self.__memmap.shape,dtype='b')
        out_idx.reshape(times,lays*vars*(rows*cols+2)+4)[:,1]=hour
        out_idx.reshape(times,lays*vars*(rows*cols+2)+4)[:,2]=date
        
        self.variables['TFLAG']=ConvertCAMxTime(self.__memmap[out_idx==date].view('>i'),self.__memmap[out_idx==hour],self.dimensions['VAR'])
        
        val_shape=out_idx.reshape(times,lays*vars*(rows*cols+2)+4)[:,4:].reshape(times,lays,vars,rows*cols+2)[:,:,:,1:-1].reshape(times,lays,vars,rows,cols)
        if self.VERSION=='<4.3':
            val_shape[:,:,0,:,:]=cloud
            val_shape[:,:,1,:,:]=rain
            val_shape[:,:,2,:,:]=cod
            self.__set_var('CLOUD',out_idx==cloud)
            self.__set_var('PRECIP',out_idx==rain)
            self.__set_var('COD',out_idx==cod)
        else:
            val_shape[:,:,0,:,:]=cloud
            val_shape[:,:,1,:,:]=rain
            val_shape[:,:,2,:,:]=snow
            val_shape[:,:,3,:,:]=graupel
            val_shape[:,:,4,:,:]=cod
            self.__set_var('CLOUD',out_idx==cloud)
            self.__set_var('RAIN',out_idx==rain)
            self.__set_var('SNOW',out_idx==snow)
            self.__set_var('GRAUPEL',out_idx==graupel)
            self.__set_var('COD',out_idx==cod)
        
        buf=self.__memmap[out_idx==0].reshape(vars*times*lays+times,2)
        if not (buf[:,0]==buf[:,1]).all():
            raise ValueError,"Buffer"
        
        return self.variables[key]

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
               
    def testCR(self):
        crfile=cloud_rain('../../../../testdata/met/camx_cr.20000825.hgbpa_04km.TCEQuh1_eta.v43',65,83)
        crfile.variables['TFLAG']
        self.assert_((crfile.variables['CLOUD'].mean(0).mean(1).mean(1)==array([  1.19793136e-02,   0.00000000e+00,   0.00000000e+00,
         7.77061632e-06,   6.17987826e-05,   1.97532892e-04,
         6.59637444e-04,   4.47470834e-03,   1.83031596e-02,
         5.06042875e-02,   6.82687238e-02,   5.76790497e-02,
         3.78983393e-02,   2.46342644e-02,   1.78406648e-02,
         7.92405475e-03,   4.40968759e-03,   2.93084350e-03,
         1.52575108e-03,   9.15497891e-04,   4.99821093e-04,
         2.04717653e-04,   1.37511626e-04,   2.84616799e-05,
         1.30704029e-05,   1.00149791e-05,   1.69081777e-05,
         0.00000000e+00],dtype='f')).all())
       
if __name__ == '__main__':
    unittest.main()
