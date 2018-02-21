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
listnan=struct.unpack('>f',b'\xff\xc0\x00\x00')[0]
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
                raise ValueError("The product of cols (%d) and rows (%d) must equal cells (%d)" %  (cols,rows,self.cell_count))

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
        for k,v in decor(key).items():
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
            raise KeyError("Vertical Diffusivity file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time))
        if chkvar and k<1 or k>self.nlayers:
            raise KeyError("Vertical Diffusivity file include layers 1 thru %i; you requested %i" % (self.nlayers,k))
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

    def values(self):
        for d,t,k in self.__iter__():
            yield self.seekandread(d,t,k)
            
    def items(self):
        for d,t,k in self.__iter__():
            yield d,t,k,self.seekandread(d,t,k)
        
    def keys(self):
        for d,t in self.timerange():
            for k in range(1,self.nlayers+1):
                yield d,t,k
    __iter__=keys

    def getArray(self):
        a=zeros((self.time_step_count,self.nlayers,len(self.dimensions['ROW']),len(self.dimensions['COL'])),'f')
        for ti,(d,t) in enumerate(self.timerange()):
            for ki,k in enumerate(range(1,self.nlayers+1)):
                self.seekandreadinto(a[ti,ki,...],d,t,k)
        return a

    def timerange(self):
        return timerange((self.start_date,self.start_time),(self.end_date,self.end_time+self.time_step),self.time_step)

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testKV(self):
        import PseudoNetCDF.testcase
        vdfile=one3d(PseudoNetCDF.testcase.camxfiles_paths['vertical_diffusivity'],4,5)
        self.assert_((vdfile.variables['UNKNOWN']==array([1.00000000e+00, 4.76359320e+00, 1.92715893e+01, 1.52158489e+01, 7.20601225e+00, 1.84097159e+00, 2.63084507e+01, 1.27621298e+01, 1.06348248e+01, 2.22587357e+01, 1.00000000e+00, 1.69009724e+01, 1.00000000e+00, 2.23075104e+01, 1.27485418e+01, 1.45508013e+01, 1.45637455e+01, 2.95294094e+01, 3.02676849e+01, 2.84974957e+01, 1.00000000e+00, 6.93706131e+00, 3.07418957e+01, 3.41300621e+01, 1.66266994e+01, 1.00000000e+00, 2.53920174e+01, 1.92539787e+01, 2.32906532e+01, 5.96702042e+01, 1.00000000e+00, 2.24458847e+01, 1.00000000e+00, 5.45038452e+01, 3.45825729e+01, 3.43578224e+01, 3.76071548e+01, 6.76799850e+01, 7.33648529e+01, 7.35801239e+01, 1.00000000e+00, 6.38178444e+00, 4.39327278e+01, 5.00754166e+01, 2.44474106e+01, 1.00000000e+00, 3.51700935e+01, 2.71137428e+01, 3.40312347e+01, 8.85775909e+01, 1.00000000e+00, 3.13994522e+01, 1.00000000e+00, 8.07169266e+01, 5.12892876e+01, 5.05734329e+01, 5.56603966e+01, 1.00394188e+02, 1.08980370e+02, 1.09251083e+02, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.64916098e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.22205174e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.05408611e+01, 1.26687222e+01, 1.21386652e+01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.93040049e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.39350021e+00, 1.02697349e+00, 1.00000000e+00, 1.00000000e+00, 1.82250175e+01, 2.90407104e+01, 2.83827496e+01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 2.02609706e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.44422662e+00, 1.02998519e+00, 1.00000000e+00, 1.00000000e+00, 2.60322971e+01, 4.26534195e+01, 4.17046585e+01],dtype='f').reshape(2,3,4,5)).all())       
if __name__ == '__main__':
    unittest.main()
