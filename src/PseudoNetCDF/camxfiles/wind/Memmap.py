__all__=['wind']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- wind Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              wind files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables, OrderedDict
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class wind(PseudoNetCDFFile):
    """
    wind provides a PseudoNetCDF interface for CAMx
    wind files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> wind_path = 'camx_wind.bin'
        >>> rows,cols = 65,83
        >>> windfile = wind(wind_path,rows,cols)
        >>> windfile.variables.keys()
        ['TFLAG', 'U', 'V']
        >>> v = windfile.variables['V']
        >>> tflag = windfile.variables['TFLAG']
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
        >>> windfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    __units='m/s'
    def __init__(self,rffile,rows,cols):
        rf=OpenRecordFile(rffile)
        self.__time_hdr_fmts={12: "fii", 8: "fi"}[rf.record_size]
        self.__time_hdr_fmts_size=rf.record_size
        self.STIME,self.SDATE=rf.unpack("fi")
        
        rf.next()
        lays=1
        record_size=rf.record_size
        while rf.record_size==record_size:
            lays+=1
            rf.next()
        self.__dummy_length=(rf.record_size+8)/4
        lays/=2
        record=rows*cols*4+8
        total_size=self.__dummy_length
        times=0
        while total_size<rf.length:
            total_size+=record*2*lays+self.__time_hdr_fmts_size+8
            times+=1
        times-=1
        
        self.variables=OrderedDict
        del rf
        self.createDimension('TSTEP',times)
        self.createDimension('DATE-TIME',2)
        self.createDimension('LAY',lays)
        self.createDimension('ROW',rows)
        self.createDimension('COL',cols)
        self.createDimension('VAR',2)
        
        self.NVARS=len(self.dimensions['VAR'])
        self.NLAYS=len(self.dimensions['LAY'])
        self.NROWS=len(self.dimensions['ROW'])
        self.NCOLS=len(self.dimensions['COL'])
        self.FTYPE=1
        
        self.__memmap=memmap(rffile,'>f','r',offset=0)
        
        if self.__time_hdr_fmts_size==12:
            self.LSTAGGER=self.__memmap[3].view('i')
        else:
            self.LSTAGGER=nan

        self.variables=PseudoNetCDFVariables(self.__variables,['TFLAG','U','V'])
    
    def __variables(self,k):
        self.__add_variables()
        return self.variables[k]
    
    def __decorator(self,k,pncfv):
        decor={'TFLAG': {'units': 'DATE-TIME','long_name': 'TFLAG','var_desc':'Time flag'}}
        for k,v in decor.get(k,{'units': 'm/s','long_name':k,'var_desc':k}).iteritems():
            setattr(pncfv,k,v)
        return pncfv
        
    def __add_variables(self):
        tsteps=len(self.dimensions['TSTEP'])
        lays=len(self.dimensions['LAY'])
        rows=len(self.dimensions['ROW'])
        cols=len(self.dimensions['COL'])
        offset=len(self.__time_hdr_fmts)+2
        block=(rows*cols+2)*2*lays
        out_idx=zeros(self.__memmap.shape,'b')
        for t in range(tsteps):
            start=(t+1)*offset+t*block+t*self.__dummy_length
            stop=start+block
            
            out_idx[start:stop].reshape(lays*2,rows*cols+2)[:,1:-1]=1
            out_idx[start:stop].reshape(lays*2,rows*cols+2)[:,[0,-1]]=2
            out_idx[start-offset:start]=3

        buffer=self.__memmap[out_idx==2].reshape(tsteps,lays,2,2)
        if not (buffer[:,:,:,0]==buffer[:,:,:,1]).all():
            raise ValueError,'Fortran unformatted record start and end padding do not match.'
        date=self.__memmap[out_idx==3].reshape(tsteps,(out_idx==3).sum()/tsteps)[:,2].view('>i')
        time=self.__memmap[out_idx==3].reshape(tsteps,(out_idx==3).sum()/tsteps)[:,1]

        self.variables['TFLAG']=ConvertCAMxTime(date,time,2)
        self.variables['U']=self.__decorator('U',PseudoNetCDFVariable(self,'U','f',('TSTEP','LAY','ROW','COL'),values=self.__memmap[out_idx==1].reshape(tsteps,lays,2,rows,cols)[:,:,0,:,:]))
        self.variables['V']=self.__decorator('V',PseudoNetCDFVariable(self,'V','f',('TSTEP','LAY','ROW','COL'),values=self.__memmap[out_idx==1].reshape(tsteps,lays,2,rows,cols)[:,:,1,:,:]))
        
    
class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testWD(self):
        import PseudoNetCDF.testcase
        wdfile=wind(PseudoNetCDF.testcase.camxfiles_paths['wind'],4,5)
        wdfile.variables['TFLAG']
        self.assert_((wdfile.variables['V'][:]==array([   -1.73236704e+00, -1.99612117e+00, -3.00912833e+00, -3.92667437e+00, -3.49521232e+00, 2.04542422e+00, 8.57666790e-01, -1.71201074e+00, -4.24386787e+00, -5.37704515e+00, 1.85697508e+00, 6.34313405e-01, -1.21529281e+00, -3.03180861e+00, -4.36278439e+00, -1.90753967e-01, -1.08261776e+00, -1.73634803e+00, -2.10829663e+00, -2.28424144e+00, -1.88443780e+00, -2.02582169e+00, -3.09955978e+00, -4.14587784e+00, -3.72402787e+00, 2.16277528e+00, 8.94082963e-01, -1.86343944e+00, -4.58147812e+00, -5.81837606e+00, 1.97949493e+00, 6.12511635e-01, -1.35096896e+00, -3.25313163e+00, -4.67790413e+00, -1.89851984e-01, -1.16381800e+00, -1.84269297e+00, -2.21348834e+00, -2.40952253e+00, -2.04972148e+00, -2.11795568e+00, -3.06094027e+00, -4.11207581e+00, -3.74964952e+00, 2.09780049e+00, 8.01259458e-01, -1.90404522e+00, -4.59170580e+00, -5.83114100e+00, 1.97475648e+00, 5.54396451e-01, -1.41695607e+00, -3.28227353e+00, -4.67724609e+00, -1.94723800e-01, -1.18353117e+00, -1.86556363e+00, -2.22842574e+00, -2.42080784e+00, -1.65720737e+00, -1.58054411e+00, -2.25336742e+00, -3.06462526e+00, -2.47261453e+00, 1.37642264e+00, 1.16142654e+00, -6.82058990e-01, -2.68112469e+00, -3.38680530e+00, 1.80796599e+00, 1.48641026e+00, -1.71508826e-02, -1.68607295e+00, -2.89399385e+00, 3.40398103e-01, 3.25049832e-02, -5.91206312e-01, -1.19038010e+00, -1.52301860e+00, -1.83006203e+00, -1.74505961e+00, -2.50190806e+00, -3.29507184e+00, -2.65367699e+00, 1.55719578e+00, 1.25234461e+00, -9.14537191e-01, -3.16307521e+00, -4.00584650e+00, 2.07018161e+00, 1.60957754e+00, -1.46312386e-01, -2.04018188e+00, -3.40377665e+00, 4.11731720e-01, -2.29119677e-02, -7.27373540e-01, -1.35116744e+00, -1.70711970e+00, -1.72859466e+00, -1.73683071e+00, -2.65253377e+00, -3.43689489e+00, -2.75470304e+00, 1.58366191e+00, 1.19656324e+00, -1.09935236e+00, -3.38544369e+00, -4.26436615e+00, 2.04826832e+00, 1.53576791e+00, -2.60809243e-01, -2.18679833e+00, -3.59082842e+00, 3.77060443e-01, -1.05680525e-01, -8.10511589e-01, -1.40993130e+00, -1.76300752e+00],dtype='f').reshape(2,3,4,5)).all())
       
TestSuite=unittest.makeSuite(TestMemmap,'test')               
if __name__ == '__main__':
    unittest.main()
