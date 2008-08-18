HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['one3d']
#Distribution packages
import unittest
import struct

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

class one3d(PseudoNetCDFFile):
    """
    vertical_diffusivity is intended to be an interface to the kv
    CAMx input file
    """
    
    id_fmt="fi"
    data_fmt="f"
    var_name="UNKNOWN"
    units="UNKNOWN"
    def __init__(self,rf,rows,cols):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        
        self.rffile=rf

        self.__memmap=memmap(self.rffile,'>f','r',offset=0)

        self.__record_items=rows*cols+4

        self.__records=self.__memmap.shape[0]/self.__record_items
        time_date=array(self.__memmap.reshape(self.__records,self.__record_items)[:,1:3])

        lays=where(time_date!=time_date[newaxis,0])[0][0]

        new_hour=slice(0,None,lays)

        dates=time_date[:,1].view('>i')
        times=time_date[:,0]

        self.__tflag=array([dates[new_hour],times[new_hour]],dtype='>f').swapaxes(0,1)
        time_steps=self.__records/lays

        self.dimensions={'VAR': 1, 'TSTEP': time_steps,'COL': cols,'ROW': rows,'LAY': lays,'DATE-TIME':2}
        
        self.FTYPE=1
        self.NVARS=1
        self.NCOLS=cols
        self.NROWS=rows
        self.NLAYS=lays
        self.NTHIK=1
        
        self.variables=PseudoNetCDFVariables(self.__variables,[self.var_name])
        v=self.variables['TFLAG']=ConvertCAMxTime(self.__tflag[:,0],self.__tflag[:,1],1)
        self.SDATE,self.STIME=v[0,0,:]

    def __decorator(self,name,pncfv):
        decor=lambda *args: dict(units=self.units, var_desc=self.var_name.ljust(16), long_name=self.var_name.ljust(16))
        for k,v in decor(name).iteritems():
            setattr(pncfv,k,v)
        return pncfv
        
    def __variables(self,k):
        tsteps=self.dimensions['TSTEP']
        lays=self.dimensions['LAY']
        rows=self.dimensions['ROW']
        cols=self.dimensions['COL']
        return self.__decorator(k,PseudoNetCDFVariable(self,k,'f',('TSTEP','LAY','ROW','COL'),self.__memmap.reshape(self.__records,self.__record_items)[:,3:-1].reshape(tsteps,lays,rows,cols)))

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testKV(self):
        vdfile=one3d('../../../../testdata/met/camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke',65,83)
        vdfile.variables['TFLAG']
        self.assert_((vdfile.variables['UNKNOWN'].mean(0).mean(1).mean(1)==array([  13.65080357,   34.39198303,   68.02783966,   95.5898819 ,
          109.25765991,  112.92014313,  108.32209778,   97.25794983,
          84.1328125 ,   65.92033386,   46.97774506,   25.8343792 ,
          9.80327034,    2.89653206,    1.26993668,    1.12098336,
          1.13557184,    1.13372564,    1.19559622,    1.1675849 ,
          1.18877947,    1.18713808,    1.02371764,    1.02544105,
          1.21638143,    1.34624374,    1.03213251,    1.        ],dtype='f')).all())
    
       
if __name__ == '__main__':
    unittest.main()
