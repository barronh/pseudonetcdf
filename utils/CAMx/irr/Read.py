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
try:
    from Scientific.IO.NetCDF import NetCDFFile as ncf
except:
    from pynetcdf import NetCDFFile as ncf

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd,timerange
from pyPA.utils.util import cartesian,sliceit
from pyPA.utils.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables


#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

__all__=['irr']

class irr(PseudoNetCDFFile):
    """
    irr is intended to be an interface to the integrated
    reaction rate file produced from CAMx with PA turned on
    """
    id_fmt="if5i"
    data_fmt="f"
    def __init__(self,rf,units='umol/hr',conva=None):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.rffile=OpenRecordFile(rf)
        self.rffile.infile.seek(0,2)
        if self.rffile.infile.tell()<2147483648L:
            warn("For greater speed on files <2GB use ipr_memmap")
        self.rffile.infile.seek(0,0)
        self.__readheader()
        self.__gettimestep()
        self.units=units
        #__conv is a conversion array that comes from ipr
        #if units is not umol/hr, conv must be provided
        self.__conv=conva
        if self.units!='umol/hr' and self.__conv==None:
            raise ValueError, "When units are provided, a conversion array dim(t,z,x,y) must also be provided"
        varkeys=['IRR_%d' % i for i in range(1,self.nrxns+1)]

        domain=self.padomains[0]
        self.dimensions=dict(TSTEP=self.time_step_count,COL=domain['iend']-domain['istart']+1,ROW=domain['jend']-domain['jstart']+1,LAY=domain['tlay']-domain['blay']+1)
        self.variables=PseudoNetCDFVariables(self.__var_get,varkeys)

    def __var_get(self,key):
        constr=lambda k: self.__variables(k,pagrid=0)
        decor=lambda k: dict(units='ppm/hr', var_desc=k.ljust(16), long_name=k.ljust(16))
        values=constr(key)
        
        var=self.createVariable(key,'f',('TSTEP','LAY','ROW','COL'))
        var.assignValue(values)
        for k,v in decor(key).iteritems():
            setattr(var,k,v)
        return var

    def __variables(self,rxn,pagrid=0):
        rxn=int(rxn[4:])
        return self.getArray(pagrid=pagrid,nrxns=rxn-1).squeeze()

    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        self.runmessage=self.rffile.read("80s")
        self.start_date,self.start_time,self.end_date,self.end_time=self.rffile.read("ifif")
        self.SDATE=self.start_date
        self.STIME=self.start_time
        
        self.grids=[]
        for grid in range(self.rffile.read("i")[-1]):
            self.grids.append(
                            dict(
                                zip(
                                    ['orgx','orgy','ncol','nrow','xsize','ysize','iutm'], 
                                    self.rffile.read("iiiiiii")
                                    )
                                )
                            )
        
        self.padomains=[]
        for padomain in range(self.rffile.read("i")[-1]):
            self.padomains.append(
                                dict(
                                    zip(
                                        ['grid','istart','iend','jstart','jend','blay','tlay'],
                                        self.rffile.read("iiiiiii")
                                        )
                                    )
                                )
        self.nrxns=self.rffile.read('i')[-1]
        
        self.data_start_byte=self.rffile.record_start
        self.record_fmt=self.id_fmt + str(self.nrxns) + self.data_fmt
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        
    def __gettimestep(self):
        """
        Header information provides start and end date, but does not
        indicate the increment between.  This routine reads the first
        and second date/time and initializes variables indicating the
        timestep length and the anticipated number.
        """
        self.activedomain=self.padomains[0]
        self.rffile._newrecord(
                        self.data_start_byte+(
                                    self.__jrecords(0,self.padomains[0]['jend'])*
                                    self.padded_size
                                    )
                        )
        date,time=self.rffile.read("if")
        self.time_step=timediff((self.start_date,self.start_time),(date,time))
        self.TSTEP=self.time_step
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time))/self.time_step)

    def __gridrecords(self,pagrid):
        """
        routine returns the number of records to increment from the
        data start byte to find the pagrid
        """
        ntime=self.__timerecords(pagrid,(self.end_date,int(self.end_time+self.time_step)))
        return ntime
        
    def __timerecords(self,pagrid,(d,t)):
        """
        routine returns the number of records to increment from the
        data start byte to find the first time
        """
        nsteps=int(timediff((self.start_date,self.start_time+self.time_step),(d,t))/self.time_step)
        nj=self.__jrecords(pagrid,self.padomains[pagrid]['jend']+1)
        return nsteps*nj
        
    def __irecords(self,pagrid,i):
        """
        routine returns the number of records to increment from the
        data start byte to find the first icell
        """
        ni=i-self.activedomain['istart']
        nk=self.__krecords(pagrid,self.activedomain['tlay']+1)
        return ni*nk
        
    def __jrecords(self,pagrid,j):
        """
        routine returns the number of records to increment from the
        data start byte to find the first jcell
        """
        nj=j-self.activedomain['jstart']
        ni=self.__irecords(pagrid,self.activedomain['iend']+1)
        return nj*ni
        
    def __krecords(self,pagrid,k):
        """
        routine returns the number of records to increment from the
        data start byte to find the first kcell
        """
        return k-self.activedomain['blay']

    def __recordposition(self,pagrid,date,time,i,j,k):
        """ 
        routine uses pagridrecords, timerecords,irecords,
        jrecords, and krecords multiplied by the fortran padded size
        to return the byte position of the specified record
        
        pagrid - integer
        date - integer
        time - float
        i - integer
        j - integer
        k - integer
        """
        records=0
        for pag in range(pagrid):
            records+=__gridrecords(pag)
        records+=self.__timerecords(pagrid,(date,time))
        records+=self.__jrecords(pagrid,j)
        records+=self.__irecords(pagrid,i)
        records+=self.__krecords(pagrid,k)
        return records*self.padded_size+self.data_start_byte
        
    def seek(self,pagrid=1,date=None,time=None,i=1,j=1,k=1):
        """
        Move file cursor to beginning of specified record
        see __recordposition for a definition of variables
        """
        if date==None:
            date=self.start_date
        if time==None:
            time=self.start_time
        self.activedomain=self.padomains[pagrid]
        #if pagrid>=len(self.padomains):
        #    raise KeyError, "IRR file contains %i PA domains; you requested the %i" % (len(self.padomains),pagrid+1)
        #if timediff((self.end_date,self.end_time),(date,time))>0 or timediff((self.start_date,self.start_time),(date,time))<0:
        #    raise KeyError, "IRR file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time)
        #if i<self.activedomain['istart'] or i>self.activedomain['iend']:
        #    raise KeyError, "IRR file i indexes include %i thru %i; you requested %i" % (self.activedomain['istart'],self.activedomain['iend'],i)
        #if j<self.activedomain['jstart'] or j>self.activedomain['jend']:
        #    raise KeyError, "IRR file j indexes include %i thru %i; you requested %i" % (self.activedomain['jstart'],self.activedomain['jend'],j)
        #if k<self.activedomain['blay'] or k>self.activedomain['tlay']:
        #    raise KeyError, "IRR file k indexes include %i thru %i; you requested %i" % (self.activedomain['blay'],self.activedomain['tlay'],k)
        
        self.rffile._newrecord(self.__recordposition(pagrid,date,time,i,j,k))
    
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
    
    def seekandreadinto(self,dest,pagrid=1,date=None,time=None,i=1,j=1,k=1):
        """
        see seek and read_into
        """
        self.seek(pagrid,date,time,i,j,k)
        return self.read_into(dest)
    
    def seekandread(self,pagrid=1,date=None,time=None,i=1,j=1,k=1):
        """
        see seek and read
        """
        self.seek(pagrid,date,time,i,j,k)
        return self.read()

    def iteritems(self,pagrid=0):
        for pagrid,d,t,i,j,k in self.iterkeys(pagrid):
            return pagrid,d,t,i,j,k,self.seekandread(pagrid,d,t,i,j,k)
 
    def itervalues(self,pagrid=0):
        for pagrid,d,t,i,j,k in self.iterkeys(pagrid):
            return self.seekandread(pagrid,d,t,i,j,k)
    
    def iterkeys(self,pagrid=0):
        domain=self.padomains[pagrid]
        for d,t in self.timerange():
            for i in range(domain['istart'],domain['iend']):
                for j in range(domain['jstart'],domain['jend']):
                    for k in range(domain['kstart'],domain['kend']):
                        return pagrid,d,t,i,j,k
                        
    def getArray(self,pagrid=1,nrxns=slice(1,None),krange=slice(1,None),nx=slice(None),ny=slice(None)):
        domain=self.padomains[pagrid]
        istart=domain['istart']
        iend=domain['iend']
        jstart=domain['jstart']
        jend=domain['jend']
        kend=domain['tlay']
        
        krange=sliceit(krange)
        nrxns=sliceit(nrxns)
        nx=sliceit(nx)
        ny=sliceit(ny)
        
        a=zeros(
            (
                len([t for t in self.timerange()]),
                self.nrxns,
                len(xrange(*krange.indices(kend+1))),
                len(range(jstart,jend+1)),
                len(range(istart,iend+1))
            ),'f')
        for ti,(d,t) in enumerate(self.timerange()):
            for ii,i in enumerate(xrange(istart,iend+1)):
                for ji,j in enumerate(xrange(jstart,jend+1)):
                    for ki,k in enumerate(xrange(*krange.indices(kend+1))):
                            self.seekandreadinto(a[ti,...,ki,ji,ii],pagrid,d,t,i,j,k)
        return a[:,nrxns,:,:,:]
            
    def timerange(self):
        return timerange((self.start_date,self.start_time+self.time_step),timeadd((self.end_date,self.end_time),(0,self.time_step)),self.time_step)

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testIRR(self):
        emissfile=irr('../../../../testdata/ei/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
        self.assert_(1==2)
       
if __name__ == '__main__':
    unittest.main()
