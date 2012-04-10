__all__=['wind']
__doc__ = """
.. _Read
:mod:`Read` -- wind Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              wind files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
    
    time_hdr_fmts={12: "fii", 8: "fi"}
    data_fmt="f"
    def __init__(self,rf,rows=None,cols=None):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.rffile=OpenRecordFile(rf)
        self.time_hdr_fmt=self.time_hdr_fmts[self.rffile.record_size]
        
        self.time_hdr_size=struct.calcsize(self.time_hdr_fmt)
        self.padded_time_hdr_size=struct.calcsize("ii"+self.time_hdr_fmt)
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

        self.variables=PseudoNetCDFVariables(self.__var_get,['U','V'])

    def __var_get(self,key):
        constr=lambda uv: self.getArray()[:,{'U': 0, 'V': 1}[uv],:,:,:].copy()
        decor=lambda uv: dict(units='m/s', var_desc=uv.ljust(16), long_name=uv.ljust(16))
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
        if self.time_hdr_fmt=='fii':
            self.start_time,self.start_date,self.lstagger=self.rffile.read(self.time_hdr_fmt)
        elif self.time_hdr_fmt=='fi':
            self.start_time,self.start_date=self.rffile.read(self.time_hdr_fmt)
            self.lstagger=None
        else:
            raise NotImplementedError, "Header format is unknown"
            
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        self.cell_count=self.record_size/struct.calcsize(self.data_fmt)
        self.record_fmt=self.data_fmt*self.cell_count
        
    def __gettimestep(self):
        """
        Header information provides start and end date, but does not
        indicate the increment between.  This routine reads the first
        and second date/time and initializes variables indicating the
        timestep length and the anticipated number.
        """
        #This is a bit of a hack, but should work:
        #Search for the next record that is the same
        #length as self.padded_time_hdr_size
        #
        #This should be the next date record
        #1) date - startdate = timestep
        #2) (record_start - self.padded_time_hdr_size)/self.padded_size = klayers
        self.rffile.next()
        while not self.rffile.record_size==self.time_hdr_size:
            self.rffile.next()
        
        dist_btwn_dates=self.rffile.record_start - self.padded_time_hdr_size
        self.nlayers=(dist_btwn_dates)/self.padded_size/2
        
        if self.time_hdr_fmt=="fi":
            time,date=self.rffile.read(self.time_hdr_fmt)
        elif self.time_hdr_fmt=="fii":
            time,date,lstagger=self.rffile.read(self.time_hdr_fmt)

        self.time_step=timediff((self.start_date,self.start_time),(date,time))
        
        while True:
            try:
                self.rffile._newrecord(self.rffile.record_start+dist_btwn_dates)
                self.rffile.tell()
                if self.time_hdr_fmt=="fi":
                    self.end_time,self.end_date=self.rffile.read(self.time_hdr_fmt)
                elif self.time_hdr_fmt=="fii":
                    self.end_time,self.end_date,lstagger=self.rffile.read(self.time_hdr_fmt)
            except:
                break
        
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time))/self.time_step)+1
        
    def __layerrecords(self,k):
        return k-1
        
    def __timerecords(self,(d,t)):
        """
        routine returns the number of records to increment from the
        data start byte to find the first time
        """
        nsteps=int(timediff((self.start_date,self.start_time),(d,t))/self.time_step)
        nlays=self.__layerrecords(self.nlayers+1)
        return nsteps*nlays
        
    def __recordposition(self,date,time,k,duv):
        """ 
        routine uses pagridrecords, timerecords,irecords,
        jrecords, and krecords multiplied by the fortran padded size
        to return the byte position of the specified record
        
        pagrid - integer
        date - integer
        time - float
        duv - integer (0=date,1=uwind,2=vwind)
        """
        bytes=self.data_start_byte
        nsteps=self.__timerecords((date,time))
        bytes+=int(nsteps/self.nlayers)*self.padded_time_hdr_size
        bytes+=int(nsteps/self.nlayers)*12
        bytes+=nsteps*self.padded_size*2
        if not duv==0:
            bytes+=self.padded_time_hdr_size
            bytes+=self.__layerrecords(k)*2*self.padded_size
        if duv==2:
            bytes+=self.padded_size
        return bytes
        
    def seek(self,date=None,time=None,k=1,uv=1):
        """
        Move file cursor to beginning of specified record
        see __recordposition for a definition of variables
        """
        if date==None:
            date=self.start_date
        if time==None:
            time=self.start_time
        chkvar=True
        if chkvar and timediff((self.end_date,self.end_time),(date,time))>0 or timediff((self.start_date,self.start_time),(date,time))<0:
            raise KeyError, "Wind file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time)
        if chkvar and uv<0 or uv >2:
            raise KeyError, "Wind file includes Date (uv: 0), u velocity (uv: 1) and v velocity (uv: 2); you requested %i" % (uv)
        
        self.rffile._newrecord(self.__recordposition(date,time,k,uv))
        
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
        return read_into(self.rffile,dest,"",self.data_fmt)
        
    def seekandreadinto(self,dest,date=None,time=None,k=1,duv=1):
        """
        see seek and read_into
        """
        self.seek(date,time,k,duv)
        self.read_into(dest)
        
    def seekandread(self,date=None,time=None,k=1,duv=1):
        """
        see seek and read
        """
        self.seek(date,time,k,duv)
        return self.read()
        
    def iterkeys(self):
        for d,t in timerange((self.start_date,self.start_time),timeadd((self.end_date,self.end_time),(0,self.time_step)),self.time_step):
            for k in range(1,self.nlayers+1):
                yield d,t,k
                
    def itervalues(self):
        for d,t,k in self.iterkeys():
            yield self.seekandread(d,t,k,1),self.seekandread(d,t,k,2)
    
    def iteritems(self):
        for d,t,k in self.iterkeys():
            yield d,t,k,self.seekandread(d,t,k,1),self.seekandread(d,t,k,2)
    
    __iter__=iterkeys
    
    def getArray(self,krange=slice(1,None)):
        if type(krange) != slice :
            if type(krange)==tuple:
                krange = slice(*krange)
            if type(krange)==int:
                krange=slice(krange,krange+1)
        a=zeros(
            (
                self.time_step_count ,
                2 ,
                len(xrange(*krange.indices(self.nlayers+1))),
                len(self.dimensions['ROW']),
                len(self.dimensions['COL']),
            ),'f')
        for i,(d,t) in enumerate(self.timerange()):
            for uv in range(1,3):
                for ki,k in enumerate(xrange(*krange.indices(self.nlayers+1))):
                    self.seekandreadinto(a[i,uv-1,k-1,:,:],d,t,k,uv)
        return a
    
    def timerange(self):
        return timerange((self.start_date,self.start_time),timeadd((self.end_date,self.end_time),(0,self.time_step)),self.time_step)


class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testWD(self):
        wdfile=wind('../../../../testdata/met/camx_wind.20000825.hgbpa_04km.TCEQuh1_eta.v43',65,83)
        self.assert_((wdfile.variables['V'].mean(0).mean(1).mean(1)==array([ 1.11006343,  1.55036175,  2.00059319,  2.1433928 ,  2.22216082,
        2.32335973,  2.35524583,  2.35921693,  2.33058643,  2.28071952,
        2.23793173,  2.04576039,  2.01241875,  2.30642509,  2.36914039,
        2.21509433,  1.73041999,  1.55570471,  1.25667989,  1.12686849,
        0.8800422 ,  1.38759625,  1.08145201, -0.32430261, -2.38370752,
       -4.76881075, -7.09392786, -5.76552343],dtype='f')).all())
       
        
if __name__ == '__main__':
    unittest.main()
