__all__=['temperature']
__doc__ = """
.. _Read
:mod:`Read` -- temperature Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              temperature files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
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

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables


#for use in identifying uncaught nan
listnan=struct.unpack('>f',b'\xff\xc0\x00\x00')[0]
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
        self.rffile=OpenRecordFile(rf)
        self.id_size=struct.calcsize(self.id_fmt)
        self.__readheader()
        self.__gettimestep()
        if rows==None and cols==None:
            rows=self.cell_count
            cols=1
        elif rows==None:
            rows=self.cell_count/cols
        elif cols==None:
            cols=self.cell_count/rows
        else:
            if cols*rows!=self.cell_count:
                raise ValueError("The product of cols (%d) and rows (%d) must equal cells (%d)" %  (cols,rows,self.cell_count))

        self.createDimension('TSTEP', self.time_step_count)
        self.createDimension('COL', cols)
        self.createDimension('ROW', rows)
        self.createDimension('LAY', self.nlayers)
        self.createDimension('SURF', 1)
        self.variables=PseudoNetCDFVariables(self.__var_get,['AIRTEMP','SURFTEMP'])

    def __var_get(self,key):
        decor=lambda k: dict(units='K',var_desc=k.ljust(16),long_name=k.ljust(16))
        constr=lambda k: self.__variables(k)
        values=constr(key)
        dims={'AIRTEMP':('TSTEP','LAY','ROW','COL'),'SURFTEMP':('TSTEP','SURF','ROW','COL')}[key]
        var=self.createVariable(key,'f',dims)
        var[:] = values
        for k,v in decor(key).items():
            setattr(var,k,v)
        return var

    def __readheader(self):
        self.data_start_byte=0
        self.rffile._newrecord(0)
        
        self.area_size=self.rffile.record_size
        self.area_count=(self.area_size-self.id_size)//struct.calcsize(self.data_fmt)
        self.area_padded_size=self.area_size+8
        self.area_fmt=self.id_fmt+self.data_fmt*(self.area_count)

        self.start_time,self.start_date=self.rffile.read(self.id_fmt)
        
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        self.cell_count=(self.record_size-self.id_size)//struct.calcsize(self.data_fmt)
        
        self.record_fmt=self.id_fmt+self.data_fmt*(self.cell_count)
    
    def __gettimestep(self):
        d,t=date,time=self.start_date,self.start_time
        self.nlayers=-1
        while (d,t)==(date,time):
            self.nlayers+=1
            t,d=self.rffile.read(self.id_fmt)
        self.time_step=timediff((self.start_date,self.start_time),(d,t))
        self.rffile.infile.seek(0,2)
        self.rffile.previous()
        self.end_time,self.end_date=self.rffile.read(self.id_fmt)
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time))//self.time_step)+1
    
    def __variables(self,k):
        if k=='SURFTEMP':
            out=zeros((len(self.dimensions['TSTEP']),1,len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
            vars=self.__surfmaps()
        elif k=='AIRTEMP':
            out=zeros((len(self.dimensions['TSTEP']),len(self.dimensions['LAY']),len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
            vars=self.__airmaps()
        dts = self.timerange()
        for i, v in enumerate(vars):
            out[i,...]=v
        return out
        
    def __surfpos(self):
        pos=self.data_start_byte+12
        inc=self.area_padded_size+self.padded_size*self.nlayers
        self.rffile.infile.seek(0,2)
        rflen=self.rffile.tell()
        while pos<rflen:
            yield pos
            pos+=inc
        raise StopIteration
        
    def __surfmaps(self):
        for pos in self.__surfpos():
            yield memmap(self.rffile.infile.name,'>f','r',pos,(self.area_count,)).reshape(len(self.dimensions['ROW']),len(self.dimensions['COL']))
            
    def __airpos(self):
        pos=self.area_padded_size+self.data_start_byte
        inc=self.area_padded_size+self.padded_size*self.nlayers
        self.rffile.infile.seek(0,2)
        rflen=self.rffile.tell()
        while pos<rflen:
            yield pos
            pos+=inc
        raise StopIteration
    
    def __airmaps(self):
        for pos in self.__airpos():
            yield memmap(self.rffile.infile.name,'>f','r',pos,((self.cell_count+4)*self.nlayers,)).reshape(self.nlayers,self.cell_count+4)[:,3:-1].reshape(len(self.dimensions['LAY']),len(self.dimensions['ROW']),len(self.dimensions['COL']))

    def timerange(self):
        return timerange((self.start_date,self.start_time),timeadd((self.end_date,self.end_time),(0,self.time_step),(2400,24)[int(self.time_step % 2)]),self.time_step,(2400,24)[int(self.time_step % 2)])

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
    
    def testTEMP(self):
        import PseudoNetCDF.testcase
        tempfile=temperature(PseudoNetCDF.testcase.camxfiles_paths['temperature'],4,5)
        self.assert_((tempfile.variables['AIRTEMP']==array([2.97762360e+02, 2.97261993e+02, 3.00761200e+02, 3.03811005e+02, 3.04561218e+02, 2.96350311e+02, 2.96676544e+02, 3.00992096e+02, 3.05474762e+02, 3.07840637e+02, 2.99522430e+02, 3.00271698e+02, 3.03738403e+02, 3.07201843e+02, 3.08288422e+02, 3.02957214e+02, 3.04927643e+02, 3.06630157e+02, 3.07726074e+02, 3.07380707e+02, 2.97516449e+02, 2.96920105e+02, 3.00340576e+02, 3.03413177e+02, 3.04202728e+02, 2.96074036e+02, 2.96250641e+02, 3.00632294e+02, 3.05113647e+02, 3.07390533e+02, 2.99310059e+02, 2.99901031e+02, 3.03344666e+02, 3.06782135e+02, 3.07819946e+02, 3.02657013e+02, 3.04522675e+02, 3.06167206e+02, 3.07235107e+02, 3.06883484e+02, 2.97677338e+02, 2.96919098e+02, 3.00031250e+02, 3.03082672e+02, 3.03850861e+02, 2.96460999e+02, 2.95947815e+02, 3.00303680e+02, 3.04781982e+02, 3.07048492e+02, 2.99246979e+02, 2.99508667e+02, 3.02997650e+02, 3.06450500e+02, 3.07478485e+02, 3.02246765e+02, 3.04192139e+02, 3.05832489e+02, 3.06897644e+02, 3.06546173e+02, 2.97428253e+02, 2.97174896e+02, 3.00208191e+02, 3.03096893e+02, 3.04174133e+02, 2.96558685e+02, 2.96706177e+02, 3.00862610e+02, 3.04807037e+02, 3.06937347e+02, 2.98850220e+02, 2.99482727e+02, 3.03085022e+02, 3.06456787e+02, 3.07406586e+02, 3.01888580e+02, 3.03996735e+02, 3.05916962e+02, 3.07113647e+02, 3.06539337e+02, 2.97645966e+02, 2.97326630e+02, 3.00117950e+02, 3.02804077e+02, 3.03801544e+02, 2.96783661e+02, 2.96694946e+02, 3.00722931e+02, 3.04501587e+02, 3.06560150e+02, 2.98854828e+02, 2.99314972e+02, 3.02861023e+02, 3.06150177e+02, 3.07073944e+02, 3.01700745e+02, 3.03746124e+02, 3.05626617e+02, 3.06770447e+02, 3.06172394e+02, 2.97927094e+02, 2.97691681e+02, 3.00104675e+02, 3.02464874e+02, 3.03398926e+02, 2.97336578e+02, 2.97074127e+02, 3.00716736e+02, 3.04132446e+02, 3.06129700e+02, 2.98817017e+02, 2.99221039e+02, 3.02649536e+02, 3.05787415e+02, 3.06698334e+02, 3.01333618e+02, 3.03411346e+02, 3.05317505e+02, 3.06446869e+02, 3.05815948e+02],dtype='f').reshape(2,3,4,5)).all())    
        
if __name__ == '__main__':
    unittest.main()
