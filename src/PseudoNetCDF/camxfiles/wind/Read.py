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
listnan=struct.unpack('>f',b'\xff\xc0\x00\x00')[0]
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
                raise ValueError("The product of cols (%d) and rows (%d) must equal cells (%d)" %  (cols,rows,self.cell_count))
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
        if self.time_hdr_fmt=='fii':
            self.start_time,self.start_date,self.lstagger=self.rffile.read(self.time_hdr_fmt)
        elif self.time_hdr_fmt=='fi':
            self.start_time,self.start_date=self.rffile.read(self.time_hdr_fmt)
            self.lstagger=None
        else:
            raise NotImplementedError("Header format is unknown")
            
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        self.cell_count=self.record_size//struct.calcsize(self.data_fmt)
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
        self.rffile._newrecord(0)
        self.rffile.next()
        nlayers = 0
        while not self.rffile.record_size==self.time_hdr_size:
            self.rffile.next()
            nlayers += 1
        
        self.nlayers = (nlayers-1)//2
        
        if self.time_hdr_fmt=="fi":
            time,date=self.rffile.read(self.time_hdr_fmt)
        elif self.time_hdr_fmt=="fii":
            time,date,lstagger=self.rffile.read(self.time_hdr_fmt)
        
        self.end_time,self.end_date = time, date
        
        self.time_step=timediff((self.start_date,self.start_time),(date,time))
        
        while True:
            try:
                for i in range(self.nlayers*2 + 1):
                    self.rffile.next()
                if self.rffile.record_size == 8:
                    self.end_time,self.end_date=self.rffile.read("fi")
                elif self.rffile.record_size == 12:
                    self.end_time,self.end_date,lstagger=self.rffile.read("fii")
                else:
                    raise KeyError()
            except:
                break
        
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time))/self.time_step)+1
        
    def __layerrecords(self,k):
        return k-1
        
    def __timerecords(self,dt):
        """
        routine returns the number of records to increment from the
        data start byte to find the first time
        """
        d, t = dt
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
            raise KeyError("Wind file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time))
        if chkvar and uv<0 or uv >2:
            raise KeyError("Wind file includes Date (uv: 0), u velocity (uv: 1) and v velocity (uv: 2); you requested %i" % (uv))
        
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
        
    def keys(self):
        for d,t in timerange((self.start_date,self.start_time),timeadd((self.end_date,self.end_time),(0,self.time_step)),self.time_step):
            for k in range(1,self.nlayers+1):
                yield d,t,k
                
    def values(self):
        for d,t,k in self.keys():
            yield self.seekandread(d,t,k,1),self.seekandread(d,t,k,2)
    
    def items(self):
        for d,t,k in self.keys():
            yield d,t,k,self.seekandread(d,t,k,1),self.seekandread(d,t,k,2)
    
    __iter__=keys
    
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
                len(range(*krange.indices(self.nlayers+1))),
                len(self.dimensions['ROW']),
                len(self.dimensions['COL']),
            ),'f')
        for i,(d,t) in enumerate(self.timerange()):
            for uv in range(1,3):
                for ki,k in enumerate(range(*krange.indices(self.nlayers+1))):
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
        import PseudoNetCDF.testcase
        wdfile=wind(PseudoNetCDF.testcase.camxfiles_paths['wind'],4,5)
        self.assert_((wdfile.variables['V'][:]==array([   -1.73236704e+00, -1.99612117e+00, -3.00912833e+00, -3.92667437e+00, -3.49521232e+00, 2.04542422e+00, 8.57666790e-01, -1.71201074e+00, -4.24386787e+00, -5.37704515e+00, 1.85697508e+00, 6.34313405e-01, -1.21529281e+00, -3.03180861e+00, -4.36278439e+00, -1.90753967e-01, -1.08261776e+00, -1.73634803e+00, -2.10829663e+00, -2.28424144e+00, -1.88443780e+00, -2.02582169e+00, -3.09955978e+00, -4.14587784e+00, -3.72402787e+00, 2.16277528e+00, 8.94082963e-01, -1.86343944e+00, -4.58147812e+00, -5.81837606e+00, 1.97949493e+00, 6.12511635e-01, -1.35096896e+00, -3.25313163e+00, -4.67790413e+00, -1.89851984e-01, -1.16381800e+00, -1.84269297e+00, -2.21348834e+00, -2.40952253e+00, -2.04972148e+00, -2.11795568e+00, -3.06094027e+00, -4.11207581e+00, -3.74964952e+00, 2.09780049e+00, 8.01259458e-01, -1.90404522e+00, -4.59170580e+00, -5.83114100e+00, 1.97475648e+00, 5.54396451e-01, -1.41695607e+00, -3.28227353e+00, -4.67724609e+00, -1.94723800e-01, -1.18353117e+00, -1.86556363e+00, -2.22842574e+00, -2.42080784e+00, -1.65720737e+00, -1.58054411e+00, -2.25336742e+00, -3.06462526e+00, -2.47261453e+00, 1.37642264e+00, 1.16142654e+00, -6.82058990e-01, -2.68112469e+00, -3.38680530e+00, 1.80796599e+00, 1.48641026e+00, -1.71508826e-02, -1.68607295e+00, -2.89399385e+00, 3.40398103e-01, 3.25049832e-02, -5.91206312e-01, -1.19038010e+00, -1.52301860e+00, -1.83006203e+00, -1.74505961e+00, -2.50190806e+00, -3.29507184e+00, -2.65367699e+00, 1.55719578e+00, 1.25234461e+00, -9.14537191e-01, -3.16307521e+00, -4.00584650e+00, 2.07018161e+00, 1.60957754e+00, -1.46312386e-01, -2.04018188e+00, -3.40377665e+00, 4.11731720e-01, -2.29119677e-02, -7.27373540e-01, -1.35116744e+00, -1.70711970e+00, -1.72859466e+00, -1.73683071e+00, -2.65253377e+00, -3.43689489e+00, -2.75470304e+00, 1.58366191e+00, 1.19656324e+00, -1.09935236e+00, -3.38544369e+00, -4.26436615e+00, 2.04826832e+00, 1.53576791e+00, -2.60809243e-01, -2.18679833e+00, -3.59082842e+00, 3.77060443e-01, -1.05680525e-01, -8.10511589e-01, -1.40993130e+00, -1.76300752e+00],dtype='f').reshape(2,3,4,5)).all())       
        
if __name__ == '__main__':
    unittest.main()
