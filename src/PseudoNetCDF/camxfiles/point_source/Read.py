__all__=['point_source']
__doc__ = """
.. _Read
:mod:`Read` -- point_source Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              point_source files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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

class point_source(PseudoNetCDFFile):
    """
    point_source provides a PseudoNetCDF interface for CAMx
    point_source files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> point_source_path = 'camx_point_source.bin'
        >>> rows,cols = 65,83
        >>> point_sourcefile = point_source(point_source_path,rows,cols)
        >>> point_sourcefile.variables.keys()
        ['TFLAG', 'ETFLAG', 'TFLAG', 'XSTK', 'YSTK', 'HSTK', 'DSTK', 'TSTK',
         'VSTK', 'KCELL', 'FLOW', 'PLMHT', 'NSTKS', 'NO', 'NO2', ...]
        >>> tflag = point_sourcefile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = point_sourcefile.variables['XSTK']
        >>> v.dimensions
        ('NSTK',)
        >>> v.shape
        (38452,)
        >>> v = point_sourcefile.variables['NO2']
        >>> v.dimensions
        ('TSTEP', 'NSTK')
        >>> v.shape
        (25, 38452)
        >>> point_sourcefile.dimensions
        {'TSTEP': 25, 'NSTK': 38452}
    """
    
    emiss_hdr_fmt="10i60i3ifif"
    grid_hdr_fmt="ffiffffiiiiifff"
    cell_hdr_fmt="iiii"
    time_hdr_fmt="ifif"
    spc_fmt="10i"
    nstk_hdr_fmt="ii"
    padded_nstk_hdr_size=struct.calcsize("ii"+nstk_hdr_fmt)
    padded_time_hdr_size=struct.calcsize("ii"+time_hdr_fmt)
    stk_hdr_fmt="ffffff"
    id_fmt="i"+spc_fmt
    id_size=struct.calcsize(id_fmt)
    data_fmt="f"
    stkprops=['XSTK','YSTK','HSTK','DSTK','TSTK','VSTK']
    stktimeprops=['KCELL','FLOW','PLMHT']

    def __init__(self,rf):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.rffile=OpenRecordFile(rf)
        self.padded_time_hdr_size=struct.calcsize(self.time_hdr_fmt+"ii")
        self.__readheader()
        self.__gettimestep()
        self.__gettimeprops()
        self.createDimension('TSTEP',self.time_step_count)
        self.createDimension('STK',self.nstk)
        varkeys=['XSTK','YSTK','HSTK','DSTK','TSTK','VSTK','KCELL','FLOW','PLMHT']+[i.strip() for i in self.spcnames]
        self.variables=PseudoNetCDFVariables(self.__var_get,varkeys)
        
    def __var_get(self,key):
        constr=self.__variables
        decor=lambda *args,**kwds: {'notread': 1}

        values=constr(key)
        if key in self.stkprops:
            var=self.createVariable(key,'f',('STK',))
        else:
            var=self.createVariable(key,'f',('TSTEP','STK'))
        var[:] = values
        for k,v in decor(key).iteritems():
            setattr(var,k,v)
        return var

    def __variables(self,k):
        if k in self.stkprops:
            return array(self.stk_props)[:,self.stkprops.index(k)]
        elif k in self.stktimeprops:
            return array(self.stk_time_props)[:,:,2:][:,:,self.stktimeprops.index(k)]
        else:
            return self.getArray()[:,self.spcnames.index(k.ljust(10)),:]
            
    def header(self):
        rdum=0.
        idum=0
        ione=1
        return [
                [self.name,self.note,ione,self.nspec,self.start_date,self.start_time,self.end_date,self.end_time],
                [rdum,rdum,self.iutm,self.xorg,self.yorg,self.delx,self.dely,self.nx,self.ny,self.nz,idum,idum,rdum,rdum,rdum],
                [ione,ione,self.nx,self.ny],
                self.spcnames,
                [ione,self.nstk],
                self.stk_props,
                self.stk_time_props
                ]

    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        vals=self.rffile.read(self.emiss_hdr_fmt)
        self.name,self.note,ione,self.nspec,self.start_date,self.start_time,self.end_date,self.end_time=vals[0:10],vals[10:70],vals[70],vals[71],vals[72],vals[73],vals[74],vals[75]
        rdum,rdum,self.iutm,self.xorg,self.yorg,self.delx,self.dely,self.nx,self.ny,self.nz,idum,idum,rdum,rdum,rdum=self.rffile.read(self.grid_hdr_fmt)
        if self.nz==0:
            #Special case of gridded emissions
            #Seems to be same as avrg
            self.nlayers=1
        else:
            self.nlayers=self.nz
        ione,ione,nx,ny=self.rffile.read(self.cell_hdr_fmt)
        if not (self.nx,self.ny)==(nx,ny):
            raise ValueError, "nx, ny defined first as %i, %i and then as %i, %i" % (self.nx,self.ny,nx,ny)
        species_temp=self.rffile.read(self.nspec*self.spc_fmt)
        self.spcnames=[]
        for i in range(0,self.nspec*10,10):
            self.spcnames.append(Int2Asc(species_temp[i:i+10]))
        
        ione,self.nstk=self.rffile.read(self.nstk_hdr_fmt)

        stkprms=zeros((self.nstk*len(self.stk_hdr_fmt),),'f')
        read_into(self.rffile,stkprms,'')
        self.rffile.next()
        #self.rffile.previous()
        #self.tmplist=self.rffile.read('ffffff'*self.nstk)
        
        stkprms=stkprms.reshape((self.nstk,len(self.stk_hdr_fmt)))
        for i in range(stkprms.shape[0]):
            if stkprms[i,-1]==array_nan:
                stkprms[i,-1]=float('-nan')
        self.stk_props=stkprms.tolist()
        self.data_start_byte=self.rffile.record_start
        self.start_date,self.start_time,end_date,end_time=self.rffile.read(self.time_hdr_fmt)
        
        self.time_step=timediff((self.start_date,self.start_time),(end_date,end_time))
        self.end_time += self.time_step
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time),(2400,24)[int(self.time_step % 2)])/self.time_step)
        
        self.stk_time_prop_fmt=""+("iiiff"*self.nstk)
        self.padded_stk_time_prop_size=struct.calcsize("ii"+self.stk_time_prop_fmt)
        
        self.record_fmt=("i10i")+self.data_fmt*(self.nstk)
        self.record_size=struct.calcsize(self.record_fmt)
        self.padded_size=self.record_size+8

    def __gettimestep(self):
        """
        this is taken care of in the readheader routine
        record format provides start and end for each hour,
        which translates to t1 and t2
        """
        pass
    
    def __gettimeprops(self):
        self.stk_time_props=[]
        for ti,(d,t) in enumerate(timerange((self.start_date,self.start_time),(self.end_date,self.end_time),self.time_step,(2400,24)[int(self.time_step % 2)])):
            tmpprop=zeros((len(self.stk_time_prop_fmt)),'f')
            tmpprop[...]=self.seekandread(d,t,1,True,self.stk_time_prop_fmt)
            tmpprop=tmpprop.reshape(self.nstk,5)
            for i in range(tmpprop.shape[0]):
                if tmpprop[i,-2]==array_nan:
                    tmpprop[i,-2]=float('-nan')

            self.stk_time_props.append(tmpprop.tolist())
            
    def __timerecords(self,(d,t)):
        """
        Calculate the number of records to increment to reach time (d,t)
        """
        nsteps=timediff((self.start_date,self.start_time),(d,t),(2400,24)[int(self.time_step % 2)])
        nspec=self.__spcrecords(self.nspec+1)
        return nsteps*(nspec)
        
    def __spcrecords(self,spc):
        """
        Calculated number of records before spc
        """
        
        return spc-1

    def __recordposition(self,date,time,spc,offset=False):
        """
        Use time (d,t), spc, and k to calculate number of records before
        desired record
        
        date - integer julian
        time - float
        spc - integer
        """
        ntime=self.__timerecords((date,time))
        nhdr=((ntime/self.__spcrecords(self.nspec+1))+1)
        nspc=self.__spcrecords(spc)
        noffset=-abs(int(offset))
        byte=self.data_start_byte
        byte+=nhdr*(self.padded_time_hdr_size+self.padded_nstk_hdr_size+self.padded_stk_time_prop_size)
        byte+=(ntime+nspc)*self.padded_size
        byte+=noffset*self.padded_stk_time_prop_size
        return byte
        
    def seek(self,date=None,time=None,spc=-1,offset=False):
        """
        Move file cursor to the beginning of the specified record
        see __recordposition for parameter definitions
        """
        #chkvar=True
        #if chkvar and timediff((self.end_date,self.end_time),(date,time),24)>0 or timediff((self.start_date,self.start_time),(date,time),24)<0:
        #    raise KeyError, "Point emission file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time)
        #if chkvar and spc<1 or spc>self.nspec:
        #    raise KeyError, "Point emission file include species 1 thru %i; you requested %i" % (self.nspec,spc)

        self.rffile._newrecord(self.__recordposition(date,time,spc,offset))
    
    def read(self,fmt=None):
        """
        Provide direct access to record file read
        """
        if fmt==None:
            fmt=self.record_fmt
        return self.rffile.read(fmt)
        
    def read_into(self,dest):
        """
        Transfer values from current record to dest
        dest - numeric or numpy array
        """
        
        return read_into(self.rffile,dest,self.id_fmt,self.data_fmt)
    
    def seekandreadinto(self,dest,date=None,time=None,spc=1):
        """
        see seek and read_into
        """
        
        self.seek(date,time,spc)
        self.read_into(dest)
        
    def seekandread(self,date=None,time=None,spc=1,offset=False,fmt=None):
        """
        see seek and read
        """
        self.seek(date,time,spc,offset)
        return self.read(fmt)

    def itervalues(self):
        for d,t,spc in self.__iter__():
            yield self.seekandread(d,t,spc)
            
    def iteritems(self):
        for d,t,spc in self.__iter__():
            yield d,t,spc,self.seekandread(d,t,spc)
            
    def iterkeys(self):
        for ti,(d,t) in enumerate(self.timerange()):
            for spc in range(1,len(self.spcnames)+1):
                yield d,t,spc

    __iter__=iterkeys
    
    def getArray(self):
        a=zeros((self.time_step_count,self.nspec,self.nstk),'f')
        for ti,(d,t) in enumerate(self.timerange()):
            for spc in range(1,len(self.spcnames)+1):
                self.seekandreadinto(a[ti,spc-1,...],d,t,spc)
        return a.copy()

    def timerange(self):
        return timerange((self.start_date,self.start_time),(self.end_date,self.end_time),self.time_step,eod=24)

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testPT(self):
        emissfile=point_source('../../../../testdata/ei/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
        self.assert_(1==2)
       
if __name__ == '__main__':
    unittest.main()
