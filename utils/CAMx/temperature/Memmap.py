HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['temperature']
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

class temperature(PseudoNetCDFFile):
    id_fmt='fi'
    data_fmt='f'
    def __init__(self,rf,rows=None,cols=None):
        self.__memmap=memmap(rf,'>f','r',offset=0)
        self.dimensions={}
        rowsXcols=self.__memmap[0].view('i')/4-2
        record_length=rowsXcols+4
        records=self.__memmap.size/record_length
        
        times=self.__memmap.reshape(records,record_length)[:,1:3]
        self.STIME,self.SDATE=times[0]
        for i,(t,d) in enumerate(times):
            if (t,d)!=(self.STIME,self.SDATE):
                break
        self.SDATE=self.SDATE.view('i')
        self.createDimension('LAY',i-1)
        self.createDimension('TSTEP',times.shape[0]/i)
        
        
        if rows==None and cols==None:
            rows=rowsXcols
            cols=1
        elif rows==None:
            rows=rowsXcols/cols
        elif cols==None:
            cols=rowsXcols/rows
        else:
            if cols*rows!=rowsXcols:
                raise ValueError, "The product of cols (%d) and rows (%d) must equal cells (%d)" %  (cols,rows,rowsXcols)
        
        self.createDimension('ROW',rows)
        self.createDimension('COL',cols)
        self.createDimension('DATE-TIME',2)
        self.createDimension('VAR',2)
        
        self.NVARS=self.dimensions['VAR']
        self.NLAYS=self.dimensions['LAY']
        self.NROWS=self.dimensions['ROW']
        self.NCOLS=self.dimensions['COL']
        self.FTYPE=1
        
        self.variables=PseudoNetCDFVariables(self.__var_get,['AIRTEMP','SURFTEMP','TFLAG'])
    def __var_get(self,key):
        lays=self.dimensions['LAY']
        times=self.dimensions['TSTEP']
        rows=self.dimensions['ROW']
        cols=self.dimensions['COL']
        surf=1
        air=2
        time=3
        date=4
        out_idx=zeros(self.__memmap.shape,dtype='b').reshape(times,lays+1,rows*cols+4)
        out_idx[:,0,3:-1]=surf
        out_idx[:,1:,3:-1]=air
        out_idx[:,:,1]=time
        out_idx[:,:,2]=date
        out_idx=out_idx.ravel()
        buf=self.__memmap[out_idx==0].reshape((lays+1)*times,2)
        if not (buf[:,0]==buf[:,1]).all():
            raise ValueError,"Buffer"
        v=self.variables['SURFTEMP']=PseudoNetCDFVariable(self,'SURFTEMP','f',('TSTEP','ROW','COL'),values=self.__memmap[out_idx==1].reshape(times,rows,cols))
        v.units='K'
        v.long_name='SURFTEMP'
        v.var_desc='SURFTEMP'
        v=self.variables['AIRTEMP']=PseudoNetCDFVariable(self,'AIRTEMP','f',('TSTEP','LAY','ROW','COL'),values=self.__memmap[out_idx==2].reshape(times,lays,rows,cols))
        v.units='K'
        v.long_name='AIRTEMP'
        v.var_desc='AIRTEMP'

        date=self.__memmap[out_idx==date].view('>i')[0:None:lays+1]
        time=self.__memmap[out_idx==time].view('>i')[0:None:lays+1]
        self.variables['TFLAG']=PseudoNetCDFVariable(self,'TFLAG','f',('TSTEP','VAR','DATE-TIME'),values=ConvertCAMxTime(date,time,2))        

        return self.variables[key]
    
class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testTEMP(self):
        import pyPA.testcase
        tempfile=temperature(pyPA.testcase.CAMxTemperature,65,83)
        tempfile.variables['TFLAG']
        self.assert_((tempfile.variables['AIRTEMP'].mean(0).mean(1).mean(1)==array([300.78424072265625, 
               300.6566162109375, 300.28036499023438, 299.61105346679688, 298.90249633789062, 
               298.17489624023438, 297.45626831054688, 296.75238037109375, 296.08837890625, 
               295.47381591796875, 294.83306884765625, 294.19754028320312, 293.681396484375, 
               293.14077758789062, 292.3377685546875, 291.16769409179688, 289.93109130859375, 
               288.52767944335938, 286.92172241210938, 284.93429565429688, 282.81634521484375, 
               279.27914428710938, 274.59872436523438, 269.8642578125, 263.99749755859375, 
               252.75810241699219, 230.50096130371094, 209.34048461914062],dtype='f')).all())
               
if __name__ == '__main__':
    unittest.main()
