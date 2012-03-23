__all__=['one3d']
__doc__ = """
.. _Read
:mod:`Read` -- one3d Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              generic 1 3D variable wind files.  See 
              PseudoNetCDF.sci_var.PseudoNetCDFFile for interface details
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

class one3d(PseudoNetCDFFile):
    """
    one3d provides a PseudoNetCDF interface for CAMx
    one3d files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> one3d_path = 'camx_one3d.bin'
        >>> rows,cols = 65,83
        >>> one3dfile = one3d(one3d_path,rows,cols)
        >>> one3dfile.variables.keys()
        ['TFLAG', 'UNKNOWN']
        >>> tflag = one3dfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = one3dfile.variables['UNKNOWN']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> one3dfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    
    id_fmt="fi"
    data_fmt="f"
    var_name="UNKNOWN"
    units="UNKNOWN"
    def __init__(self,rf,rows=None,cols=None):
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
        self.createDimension('COL', cols)
        self.createDimension('ROW', rows)
        self.createDimension('LAY', self.nlayers)

        self.variables=PseudoNetCDFVariables(self.__var_get,[self.var_name])

    def __var_get(self,key):
        constr=lambda *args, **kwds: self.getArray()
        decor=lambda *args: dict(units=self.units, var_desc=self.var_name.ljust(16), long_name=self.var_name.ljust(16))
        values=constr(key)
        
        var=self.createVariable(key,'f',('TSTEP','LAY','ROW','COL'))
        var[:] = values
        for k,v in decor(key).iteritems():
            setattr(var,k,v)
        return var

    def __readheader(self):
        """
        __readheader reads the header section of the vertical diffusivity file
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
                        self.padded_size
                        )
        d,t=self.start_date,self.start_time
        self.nlayers=0
        while timediff((self.start_date,self.start_time),(d,t))==0:
            t,d=self.rffile.read(self.id_fmt)
            self.nlayers+=1
        self.time_step=timediff((self.start_date,self.start_time),(d,t))

        while True:
            try:
                self.seek(d,t,1,False)
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
        return k-1

    def __recordposition(self,date,time,k):
        """ 
        routine uses timerecords and layerrecords multiplied 
        by the fortran padded size to return the byte position of the specified record
        
        date - integer
        time - float
        k - integer
        """
        ntime=self.__timerecords((date,time))
        nk=self.__layerrecords(k)
        return (nk+ntime)*self.padded_size+self.data_start_byte
        
    def seek(self,date=None,time=None,k=1,chkvar=True):
        """
        Move file cursor to beginning of specified record
        see __recordposition for a definition of variables
        """
        if date==None:
            date=self.start_date
        if time==None:
            time=self.start_time
            
        if chkvar and timediff((self.end_date,self.end_time),(date,time))>0 or timediff((self.start_date,self.start_time),(date,time))<0:
            raise KeyError, "Vertical Diffusivity file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time)
        if chkvar and k<1 or k>self.nlayers:
            raise KeyError, "Vertical Diffusivity file include layers 1 thru %i; you requested %i" % (self.nlayers,k)
        self.rffile._newrecord(self.__recordposition(date,time,k))
        
    def read(self):
        """
        provide direct access to the underlying RecordFile read
        method
        """
        return self.rffile.read(self.record_fmt)
        
    def read_into(self,dest):
        """
        put values from rffile read into dest
        dest - numpy or numeric array
        """
        return read_into(self.rffile,dest,self.id_fmt,self.data_fmt)
        
    def seekandreadinto(self,dest,date=None,time=None,k=1):
        """
        see seek and read_into
        """
        self.seek(date,time,k)
        return self.read_into(dest)
        
    def seekandread(self,date=None,time=None,k=1):
        """
        see seek and read
        """
        self.seek(date,time,k)
        return self.read()

    def itervalues(self):
        for d,t,k in self.__iter__():
            yield self.seekandread(d,t,k)
            
    def iteritems(self):
        for d,t,k in self.__iter__():
            yield d,t,k,self.seekandread(d,t,k)
        
    def iterkeys(self):
        for d,t in self.timerange():
            for k in range(1,self.nlayers+1):
                yield d,t,k
    __iter__=iterkeys

    def getArray(self):
        a=zeros((self.time_step_count,self.nlayers,len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
        for ti,(d,t) in enumerate(self.timerange()):
            for ki,k in enumerate(range(1,self.nlayers+1)):
                self.seekandreadinto(a[ti,ki,...,...],d,t,k)
        return a

    def timerange(self):
        return timerange((self.start_date,self.start_time),(self.end_date,self.end_time+self.time_step),self.time_step)

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testKV(self):
        vdfile=one3d('../../../../testdata/met/camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke',65,83)
        self.assert_((vdfile.variables['UNKNOWN'].mean(0).mean(1).mean(1)==array([  13.65080357,   34.39198303,   68.02783966,   95.5898819 ,
          109.25765991,  112.92014313,  108.32209778,   97.25794983,
          84.1328125 ,   65.92033386,   46.97774506,   25.8343792 ,
          9.80327034,    2.89653206,    1.26993668,    1.12098336,
          1.13557184,    1.13372564,    1.19559622,    1.1675849 ,
          1.18877947,    1.18713808,    1.02371764,    1.02544105,
          1.21638143,    1.34624374,    1.03213251,    1.        ],dtype='f')).all())
       
if __name__ == '__main__':
    unittest.main()
