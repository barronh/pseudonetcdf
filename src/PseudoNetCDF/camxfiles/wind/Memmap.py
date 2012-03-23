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
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
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
        
        self.variables={}
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
        wdfile=wind(PseudoNetCDF.testcase.CAMxWind,65,83)
        wdfile.variables['TFLAG']
        self.assert_((wdfile.variables['V'].mean(0).mean(1).mean(1)==array([ 1.11006343,  1.55036175,  2.00059319,  2.1433928 ,  2.22216082,
        2.32335973,  2.35524583,  2.35921693,  2.33058643,  2.28071952,
        2.23793173,  2.04576039,  2.01241875,  2.30642509,  2.36914039,
        2.21509433,  1.73041999,  1.55570471,  1.25667989,  1.12686849,
        0.8800422 ,  1.38759625,  1.08145201, -0.32430261, -2.38370752,
       -4.76881075, -7.09392786, -5.76552343],dtype='f')).all())
       
TestSuite=unittest.makeSuite(TestMemmap,'test')               
if __name__ == '__main__':
    unittest.main()
