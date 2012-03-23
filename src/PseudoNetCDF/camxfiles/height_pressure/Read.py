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
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
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
                raise ValueError, "The product of cols (%d) and rows (%d) must equal cells (%d)" %  (cols,rows,self.cell_count)

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
        for k,v in decor(key).iteritems():
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
        self.cell_count=(self.record_size-self.id_size)/struct.calcsize(self.data_fmt)
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
        
    def __timerecords(self,(d,t)):
        """
        routine returns the number of records to increment from the
        data start byte to find the first time
        """
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
            raise KeyError, "Vertical Diffusivity file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time)
        if chkvar and k<1 or k>self.nlayers:
            raise KeyError, "Vertical Diffusivity file include layers 1 thru %i; you requested %i" % (self.nlayers,k)
        if chkvar and hp<0 or hp >1:
            raise KeyError, "Height pressure or indexed 0 and 1; you requested %i" % (hp)
            
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
    def itervalues(self):
        for d,t,k in self.__iter__():
            yield self.seekandread(d,t,k,0),self.seekandread(d,t,k,1)
            
    def iteritems(self):
        for d,t,k in self.__iter__():
            yield d,t,k,self.seekandread(d,t,k,0),self.seekandread(d,t,k,1)
        
    def iterkeys(self):
        for d,t in self.timerange():
            for k in range(1,self.nlayers+1):
                yield d,t,k
    
    __iter__=iterkeys
    
    def getArray(self,hp):
        a=zeros((self.time_step_count,len(self.dimensions['LAY']),len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
            
        for ti,(d,t) in enumerate(self.timerange()):
            for ki,k in enumerate(xrange(1,self.nlayers+1)):
                self.seekandreadinto(a[ti,ki,...,...],d,t,k,hp)
        return a
    
    def timerange(self):
        return timerange((self.start_date,self.start_time),timeadd((self.end_date,self.end_time),(0,self.time_step),(2400,24)[int(self.time_step % 2)]),self.time_step,(2400,24)[int(self.time_step % 2)])


class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testHP(self):
        hpfile=height_pressure('../../../../testdata/met/camx_zp.20000825.hgbpa_04km.TCEQuh1_eta.v43',65,83)
        self.assert_((hpfile.variables['HGHT'].mean(0).mean(1).mean(1)==array([    33.890625  ,     84.93463135,    170.56620789,    256.90753174,
          343.97067261,    431.7701416 ,    520.31616211,    609.62701416,
          699.71228027,    790.59112549,    928.42340088,   1068.12695312,
         1209.74743652,   1353.35302734,   1597.25683594,   1847.13049316,
         2103.29101562,   2366.0949707 ,   2690.75830078,   3026.26611328,
         3373.41748047,   4106.40332031,   4898.45361328,   5836.37402344,
         6961.35107422,   9166.87207031,  13057.72558594,  15177.14160156],dtype='f')).all())

        
if __name__ == '__main__':
    unittest.main()
