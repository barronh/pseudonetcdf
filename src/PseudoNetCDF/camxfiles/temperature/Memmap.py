__all__=['temperature']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- temperature Memmap interface
=============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              temperature files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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

class temperature(PseudoNetCDFFile):
    """
    temperature provides a PseudoNetCDF interface for CAMx
    temperature files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> temperature_path = 'camx_temperature.bin'
        >>> rows,cols = 65,83
        >>> temperaturefile = temperature(temperature_path,rows,cols)
        >>> temperaturefile.variables.keys()
        ['TFLAG', 'AIRTEMP', 'SURFTEMP']
        >>> tflag = temperaturefile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = temperaturefile.variables['SURFTEMP']
        >>> v.dimensions
        ('TSTEP', 'ROW', 'COL')
        >>> v.shape
        (25, 65, 83)
        >>> v = temperaturefile.variables['AIRTEMP']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> temperaturefile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    id_fmt='fi'
    data_fmt='f'
    def __init__(self,rf,rows=None,cols=None):
        self.__memmap=memmap(rf,'>f','r',offset=0)
        
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
        
        self.NVARS=len(self.dimensions['VAR'])
        self.NLAYS=len(self.dimensions['LAY'])
        self.NROWS=len(self.dimensions['ROW'])
        self.NCOLS=len(self.dimensions['COL'])
        self.FTYPE=1
        
        self.variables=PseudoNetCDFVariables(self.__var_get,['AIRTEMP','SURFTEMP','TFLAG'])
    def __var_get(self,key):
        lays=len(self.dimensions['LAY'])
        times=len(self.dimensions['TSTEP'])
        rows=len(self.dimensions['ROW'])
        cols=len(self.dimensions['COL'])
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
        time=self.__memmap[out_idx==time].view('>f')[0:None:lays+1]
        self.variables['TFLAG']=PseudoNetCDFVariable(self,'TFLAG','f',('TSTEP','VAR','DATE-TIME'),values=ConvertCAMxTime(date,time,2))        

        return self.variables[key]
    
class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testTEMP(self):
        import PseudoNetCDF.testcase
        tempfile=temperature(PseudoNetCDF.testcase.camxfiles_paths['temperature'],4,5)
        tempfile.variables['TFLAG']
        self.assert_((tempfile.variables['AIRTEMP']==array([2.97762360e+02, 2.97261993e+02, 3.00761200e+02, 3.03811005e+02, 3.04561218e+02, 2.96350311e+02, 2.96676544e+02, 3.00992096e+02, 3.05474762e+02, 3.07840637e+02, 2.99522430e+02, 3.00271698e+02, 3.03738403e+02, 3.07201843e+02, 3.08288422e+02, 3.02957214e+02, 3.04927643e+02, 3.06630157e+02, 3.07726074e+02, 3.07380707e+02, 2.97516449e+02, 2.96920105e+02, 3.00340576e+02, 3.03413177e+02, 3.04202728e+02, 2.96074036e+02, 2.96250641e+02, 3.00632294e+02, 3.05113647e+02, 3.07390533e+02, 2.99310059e+02, 2.99901031e+02, 3.03344666e+02, 3.06782135e+02, 3.07819946e+02, 3.02657013e+02, 3.04522675e+02, 3.06167206e+02, 3.07235107e+02, 3.06883484e+02, 2.97677338e+02, 2.96919098e+02, 3.00031250e+02, 3.03082672e+02, 3.03850861e+02, 2.96460999e+02, 2.95947815e+02, 3.00303680e+02, 3.04781982e+02, 3.07048492e+02, 2.99246979e+02, 2.99508667e+02, 3.02997650e+02, 3.06450500e+02, 3.07478485e+02, 3.02246765e+02, 3.04192139e+02, 3.05832489e+02, 3.06897644e+02, 3.06546173e+02, 2.97428253e+02, 2.97174896e+02, 3.00208191e+02, 3.03096893e+02, 3.04174133e+02, 2.96558685e+02, 2.96706177e+02, 3.00862610e+02, 3.04807037e+02, 3.06937347e+02, 2.98850220e+02, 2.99482727e+02, 3.03085022e+02, 3.06456787e+02, 3.07406586e+02, 3.01888580e+02, 3.03996735e+02, 3.05916962e+02, 3.07113647e+02, 3.06539337e+02, 2.97645966e+02, 2.97326630e+02, 3.00117950e+02, 3.02804077e+02, 3.03801544e+02, 2.96783661e+02, 2.96694946e+02, 3.00722931e+02, 3.04501587e+02, 3.06560150e+02, 2.98854828e+02, 2.99314972e+02, 3.02861023e+02, 3.06150177e+02, 3.07073944e+02, 3.01700745e+02, 3.03746124e+02, 3.05626617e+02, 3.06770447e+02, 3.06172394e+02, 2.97927094e+02, 2.97691681e+02, 3.00104675e+02, 3.02464874e+02, 3.03398926e+02, 2.97336578e+02, 2.97074127e+02, 3.00716736e+02, 3.04132446e+02, 3.06129700e+02, 2.98817017e+02, 2.99221039e+02, 3.02649536e+02, 3.05787415e+02, 3.06698334e+02, 3.01333618e+02, 3.03411346e+02, 3.05317505e+02, 3.06446869e+02, 3.05815948e+02],dtype='f').reshape(2,3,4,5)).all())
               
if __name__ == '__main__':
    unittest.main()
