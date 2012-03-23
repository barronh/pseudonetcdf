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
                raise ValueError, "The product of cols (%d) and rows (%d) must equal cells (%d)" %  (cols,rows,self.cell_count)

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
        for k,v in decor(key).iteritems():
            setattr(var,k,v)
        return var

    def __readheader(self):
        self.data_start_byte=0
        self.rffile._newrecord(0)
        
        self.area_size=self.rffile.record_size
        self.area_count=(self.area_size-self.id_size)/struct.calcsize(self.data_fmt)
        self.area_padded_size=self.area_size+8
        self.area_fmt=self.id_fmt+self.data_fmt*(self.area_count)

        self.start_time,self.start_date=self.rffile.read(self.id_fmt)
        
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        self.cell_count=(self.record_size-self.id_size)/struct.calcsize(self.data_fmt)
        
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
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time))/self.time_step)+1
    
    def __variables(self,k):
        if k=='SURFTEMP':
            out=zeros((len(self.dimensions['TSTEP']),1,len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
            vars=self.__surfmaps()
        elif k=='AIRTEMP':
            out=zeros((len(self.dimensions['TSTEP']),len(self.dimensions['LAY']),len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
            vars=self.__airmaps()
        for i,(d,t) in enumerate(self.timerange()):
            out[i,...]=vars.next()
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
        vdfile=temperature('../../../../testdata/met/camx_temp.20000825.hgbpa_04km.TCEQuh1_eta.v43',65,83)
        self.assert_((vdfile.variables['AIRTEMP'].mean(0).mean(1).mean(1)==array([300.78424072265625, 
               300.6566162109375, 300.28036499023438, 299.61105346679688, 298.90249633789062, 
               298.17489624023438, 297.45626831054688, 296.75238037109375, 296.08837890625, 
               295.47381591796875, 294.83306884765625, 294.19754028320312, 293.681396484375, 
               293.14077758789062, 292.3377685546875, 291.16769409179688, 289.93109130859375, 
               288.52767944335938, 286.92172241210938, 284.93429565429688, 282.81634521484375, 
               279.27914428710938, 274.59872436523438, 269.8642578125, 263.99749755859375, 
               252.75810241699219, 230.50096130371094, 209.34048461914062],dtype='f')).all())
    
        
if __name__ == '__main__':
    unittest.main()
