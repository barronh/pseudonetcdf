__all__=['height_pressure']
__doc__ = """
.. _Read
:mod:`Read` -- height_pressure Read interface
=============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              height_pressure files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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

class height_pressure(PseudoNetCDFFile):
    """
    height_pressure provides a PseudoNetCDF interface for CAMx
    height_pressure files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> height_pressure_path = 'camx_height_pressure.bin'
        >>> rows,cols = 65,83
        >>> height_pressurefile = height_pressure(height_pressure_path,rows,cols)
        >>> height_pressurefile.variables.keys()
        ['TFLAG', 'HGHT', 'PRES']
        >>> v = height_pressurefile.variables['V']
        >>> tflag = height_pressurefile.variables['TFLAG']
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
        >>> height_pressurefile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    id_fmt="fi"
    data_fmt="f"
        
    def __init__(self,rf,rows,cols):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
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
        self.createDimension('ROW', rows)
        self.createDimension('COL', cols)
        self.createDimension('LAY', self.nlayers)
        self.variables=PseudoNetCDFVariables(self.__var_get,['HGHT','PRES'])

    def __var_get(self,key):
        constr=lambda hp: self.getArray({'HGHT': 0,'PRES': 1}[hp])
        decor=lambda hp: {'HGHT': dict(units='m',var_desc='HGHT'.ljust(16),long_name='HGHT'.ljust(16)),'PRES': dict(units='hPa',var_desc='PRES'.ljust(16),long_name='PRES'.ljust(16))}[hp]
        values=constr(key)
        var=self.createVariable(key,'f',('TSTEP','LAY','ROW','COL'))
        var[:] = values
        for k,v in decor(key).items():
            setattr(var,k,v)
        return var

    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        self.data_start_byte=0
        self.start_time,self.start_date=self.rffile.read(self.id_fmt)
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        self.cell_count=(self.record_size-self.id_size)//struct.calcsize(self.data_fmt)
        self.record_fmt=self.id_fmt+self.data_fmt*(self.cell_count)
        
    def __gettimestep(self):
        """
        Header information provides start and end date, but does not
        indicate the increment between.  This routine reads the first
        and second date/time and initializes variables indicating the
        timestep length and the anticipated number.
        """
        self.rffile._newrecord(
                        self.padded_size*2
                        )
        d,t=self.start_date,self.start_time
        self.nlayers=0
        while timediff((self.start_date,self.start_time),(d,t))==0:
            t,d=self.rffile.read(self.id_fmt)
            self.rffile.next()
            self.nlayers+=1
        self.time_step=timediff((self.start_date,self.start_time),(d,t))

        while True:
            try:
                self.seek(d,t,self.nlayers,1,False)
                d,t=timeadd((d,t),(0,self.time_step))
            except:
                break
        self.end_date,self.end_time=timeadd((d,t),(0,-self.time_step))
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time))/self.time_step)+1
        
    def __timerecords(self,dt):
        """
        routine returns the number of records to increment from the
        data start byte to find the first time
        """
        d, t = dt
        nsteps=int(timediff((self.start_date,self.start_time),(d,t))/self.time_step)
        nk=self.__layerrecords(self.nlayers+1)
        return nsteps*nk
        
    def __layerrecords(self,k):
        """
        routine returns the number of records to increment from the
        data start byte to find the first klayer
        """
        return (k-1)*2

    def __recordposition(self,date,time,k,hp):
        """ 
        routine uses timerecords and layerrecords multiplied plus hp
        by the fortran padded size to return the byte position of the specified record
        
        date - integer
        time - float
        k - integer
        hp - integer (0=h,1=p)
        """
        ntime=self.__timerecords((date,time))
        nk=self.__layerrecords(k)
        return (nk+ntime+hp)*self.padded_size+self.data_start_byte
        
    def seek(self,date=None,time=None,k=1,hp=0,chkvar=True):
        """
        Move file cursor to specified record
        """
        if date==None:
            date=self.start_date
        if time==None:
            time=self.start_time
            
        if chkvar and timediff((self.end_date,self.end_time),(date,time))>0 or timediff((self.start_date,self.start_time),(date,time))<0:
            raise KeyError("Vertical Diffusivity file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time))
        if chkvar and k<1 or k>self.nlayers:
            raise KeyError("Vertical Diffusivity file include layers 1 thru %i; you requested %i" % (self.nlayers,k))
        if chkvar and hp<0 or hp >1:
            raise KeyError("Height pressure or indexed 0 and 1; you requested %i" % (hp))
            
        self.rffile._newrecord(self.__recordposition(date,time,k,hp))
        
    def read(self):
        """
        Call recordfile read method directly
        """
        return self.rffile.read(self.record_fmt)
        
    def read_into(self,dest):
        """
        put values from rffile read into dest
        dest - numpy or numeric array
        """
        return read_into(self.rffile,dest,self.id_fmt,self.data_fmt)
        
    def seekandreadinto(self,dest,date=None,time=None,k=1,hp=0):
        """
        see seek and read
        """
        self.seek(date,time,k,hp)
        return self.read_into(dest)
        
    def seekandread(self,date=None,time=None,k=1,hp=0):
        """
        see seek and read
        """
        self.seek(date,time,k,hp)
        return self.read()
    def values(self):
        for d,t,k in self.__iter__():
            yield self.seekandread(d,t,k,0),self.seekandread(d,t,k,1)
            
    def items(self):
        for d,t,k in self.__iter__():
            yield d,t,k,self.seekandread(d,t,k,0),self.seekandread(d,t,k,1)
        
    def keys(self):
        for d,t in self.timerange():
            for k in range(1,self.nlayers+1):
                yield d,t,k
    
    __iter__=keys
    
    def getArray(self,hp):
        a=zeros((self.time_step_count,len(self.dimensions['LAY']),len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
            
        for ti,(d,t) in enumerate(self.timerange()):
            for ki,k in enumerate(range(1,self.nlayers+1)):
                self.seekandreadinto(a[ti,ki,...],d,t,k,hp)
        return a
    
    def timerange(self):
        return timerange((self.start_date,self.start_time),timeadd((self.end_date,self.end_time),(0,self.time_step),(2400,24)[int(self.time_step % 2)]),self.time_step,(2400,24)[int(self.time_step % 2)])


class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testHP(self):
        import PseudoNetCDF.testcase
        hpfile=height_pressure(PseudoNetCDF.testcase.camxfiles_paths['height_pressure'],4,5)
        self.assert_((hpfile.variables['HGHT']==array([  3.38721924e+01, 3.40657959e+01, 3.41392822e+01, 3.42358398e+01, 3.42543945e+01, 3.38868408e+01, 3.40622559e+01, 3.42358398e+01, 3.44768066e+01, 3.46112061e+01, 3.37558594e+01, 3.39323730e+01, 3.42663574e+01, 3.46854248e+01, 3.48144531e+01, 3.39472656e+01, 3.41900635e+01, 3.46160889e+01, 3.48209229e+01, 3.47874756e+01, 6.78652344e+01, 6.82532959e+01, 6.84020996e+01, 6.85950928e+01, 6.86304932e+01, 6.78945312e+01, 6.82465820e+01, 6.85941162e+01, 6.90783691e+01, 6.93474121e+01, 6.76313477e+01, 6.79859619e+01, 6.86558838e+01, 6.94960938e+01, 6.97552490e+01, 6.80159912e+01, 6.85028076e+01, 6.93570557e+01, 6.97674561e+01, 6.97009277e+01, 1.01980713e+02, 1.02563843e+02, 1.02787109e+02, 1.03077759e+02, 1.03131104e+02, 1.02022949e+02, 1.02553101e+02, 1.03076904e+02, 1.03804565e+02, 1.04208984e+02, 1.01628662e+02, 1.02162842e+02, 1.03169922e+02, 1.04433838e+02, 1.04823120e+02, 1.02206909e+02, 1.02940186e+02, 1.04224609e+02, 1.04841797e+02, 1.04740723e+02, 3.38721924e+01, 3.40657959e+01, 3.41392822e+01, 3.42358398e+01, 3.42543945e+01, 3.38868408e+01, 3.40622559e+01, 3.42358398e+01, 3.44768066e+01, 3.46112061e+01, 3.37558594e+01, 3.39323730e+01, 3.42663574e+01, 3.46854248e+01, 3.48144531e+01, 3.39472656e+01, 3.41900635e+01, 3.46160889e+01, 3.48209229e+01, 3.47874756e+01, 6.78652344e+01, 6.82532959e+01, 6.84020996e+01, 6.85950928e+01, 6.86304932e+01, 6.78945312e+01, 6.82465820e+01, 6.85941162e+01, 6.90783691e+01, 6.93474121e+01, 6.76313477e+01, 6.79859619e+01, 6.86558838e+01, 6.94960938e+01, 6.97552490e+01, 6.80159912e+01, 6.85028076e+01, 6.93570557e+01, 6.97674561e+01, 6.97009277e+01, 1.01980713e+02, 1.02563843e+02, 1.02787109e+02, 1.03077759e+02, 1.03131104e+02, 1.02022949e+02, 1.02553101e+02, 1.03076904e+02, 1.03804565e+02, 1.04208984e+02, 1.01628662e+02, 1.02162842e+02, 1.03169922e+02, 1.04433838e+02, 1.04823120e+02, 1.02206909e+02, 1.02940186e+02, 1.04224609e+02, 1.04841797e+02, 1.04740723e+02],dtype='f').reshape(2, 3, 4, 5)).all())
        
if __name__ == '__main__':
    unittest.main()
