__all__=['uamiv']
__doc__ = """
.. _Read
:mod:`Read` -- uamiv Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              uamiv files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from numpy import zeros,array,where,memmap,newaxis,dtype,fromfile

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.util import sliceit
from PseudoNetCDF.camxfiles.units import get_uamiv_units
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoIOAPIVariable, PseudoNetCDFVariables


#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class uamiv(PseudoNetCDFFile):
    """
    uamiv provides a PseudoNetCDF interface for CAMx
    uamiv files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> uamiv_path = 'camx_uamiv.bin'
        >>> rows,cols = 65,83
        >>> uamivfile = uamiv(uamiv_path,rows,cols)
        >>> uamivfile.variables.keys()
        ['TFLAG', 'O3', 'NO', 'NO2', ...]
        >>> tflag = uamivfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = uamivfile.variables['O3']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> uamivfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    
    emiss_hdr_fmt="10i60i3ifif"
    grid_hdr_fmt="ffiffffiiiiifff"
    cell_hdr_fmt="iiii"
    time_hdr_fmt="ifif"
    spc_fmt="10i"
    id_fmt="i"+spc_fmt
    id_size=struct.calcsize(id_fmt)
    data_fmt="f"
    ione=1
    idum=0
    rdum=0.
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
        self.dimensions={}
        self.createDimension('LAY', self.nlayers)
        self.createDimension('COL', self.nx)
        self.createDimension('ROW', self.ny)
        self.createDimension('TSTEP', self.time_step_count)
        self.createDimension('DATE-TIME', 2)
            
        self.variables=PseudoNetCDFVariables(self.__var_get,map(str.strip,self.spcnames))

    def __var_get(self,key):
        units = get_uamiv_units(self.name, key)
        spcnames = map(str.strip, self.spcnames)
        if self.name=='EMISSIONS ':
            constr=lambda spc: self.getArray(nspec=spcnames.index(spc)).squeeze()[:,newaxis,:,:]
            decor=lambda spc: dict(units=units, var_desc=spc, long_name=spc.ljust(16))
        else:
            constr=lambda spc: self.getArray(nspec=spcnames.index(spc)).squeeze().reshape(map(len, (self.dimensions['TSTEP'],self.dimensions['LAY'],self.dimensions['ROW'],self.dimensions['COL'])))
            decor=lambda spc: dict(units=units, var_desc=spc.ljust(16), long_name=spc.ljust(16))

        values=constr(key)
        var=self.createVariable(key,'f',('TSTEP','LAY','ROW','COL'))
        var[:] = values
        for k,v in decor(key).iteritems():
            setattr(var,k,v)
        return var

    def header(self):
        rdum=self.rdum
        idum=self.idum
        ione=self.ione
        return [
                [self.name,self.note,ione,self.nspec,self.start_date,self.start_time,self.end_date,self.end_time],
                [rdum,rdum,self.iutm,self.xorg,self.yorg,self.delx,self.dely,self.nx,self.ny,self.nz,idum,idum,rdum,rdum,rdum],
                [ione,ione,self.nx,self.ny],
                self.spcnames
                ]
                
    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        vals=self.rffile.read(self.emiss_hdr_fmt)
        self.name,self.note,ione,self.nspec,self.start_date,self.start_time,self.end_date,self.end_time=vals[0:10],vals[10:70],vals[70],vals[71],vals[72],vals[73],vals[74],vals[75]
        
        self.name=Int2Asc(self.name)
        self.note=Int2Asc(self.note)
        self.rdum,rdum,self.iutm,self.xorg,self.yorg,self.delx,self.dely,self.nx,self.ny,self.nz,idum,self.idum,rdum,rdum,rdum=self.rffile.read(self.grid_hdr_fmt)

        if self.name=='EMISSIONS ':
            #Special case of gridded emissions
            #Seems to be same as avrg
            self.nlayers=1
        else:
            self.nlayers=self.nz
        self.ione,ione,nx,ny=self.rffile.read(self.cell_hdr_fmt)
        if not (self.nx,self.ny)==(nx,ny):
            raise ValueError, "nx, ny defined first as %i, %i and then as %i, %i" % (self.nx,self.ny,nx,ny)
        species_temp=self.rffile.read(self.nspec*self.spc_fmt)
        self.spcnames=[]
        for i in range(0,self.nspec*10,10):
            self.spcnames.append(Int2Asc(species_temp[i:i+10]))
        
        self.data_start_byte=self.rffile.record_start
        start_date,start_time,end_date,end_time=self.rffile.read(self.time_hdr_fmt)
        self.time_step=timediff((start_date,start_time),(end_date,end_time))
        self.time_step_count=int(timediff((self.start_date,self.start_time),(self.end_date,self.end_time),(2400,24)[int(self.time_step % 2)])//self.time_step)
        if self.name == 'AIRQUALITY':
            self.time_step_count = 1
            self.start_date = self.end_date
        self.record_size=self.rffile.record_size
        self.padded_size=self.record_size+8
        self.cell_count=(self.record_size-struct.calcsize("i10i"))/struct.calcsize(self.data_fmt)
        self.record_fmt=("i10i")+self.data_fmt*(self.cell_count)
        
    def __gettimestep(self):
        """
        this is taken care of in the readheader routine
        record format provides start and end for each hour,
        which translates to t1 and t2
        """
        pass
    
    def __timerecords(self,(d,t)):
        """
        Calculate the number of records to increment to reach time (d,t)
        """
        nsteps=int(timediff((self.start_date,self.start_time),(d,t))/self.time_step)
        nspec=self.__spcrecords(self.nspec+1)
        return nsteps*nspec
    def __layerrecords(self,k):
        """Calculate the number of records to increment to reach layer k
        """
        return k-1
    def __spcrecords(self,spc):
        """
        Calculated number of records before spc
        """
        return (spc-1)*self.__layerrecords(self.nlayers+1)

    def __recordposition(self,date,time,spc,k):
        """
        Use time (d,t), spc, and k to calculate number of records before
        desired record
        
        date - integer julian
        time - float
        spc - integer
        k - integer
        """
        ntime=self.__timerecords((date,time))
        nk=self.__layerrecords(k)
        nid=ntime/self.nspec/self.nlayers
        nspec=0
        if spc!=0:
            nid+=1
            nspec=self.__spcrecords(spc)

        return self.data_start_byte+(nspec+nk+ntime)*self.padded_size+nid*self.padded_time_hdr_size
        
    def seek(self,date=None,time=None,spc=-1,k=0,chkvar=True):
        """
        Move file cursor to the beginning of the specified record
        see __recordposition for parameter definitions
        """
        spc+=1
        if date==None:
            date=self.start_date
        if time==None:
            time=self.start_time
            
        if chkvar and timediff((self.end_date,self.end_time),(date,time),24)>0 or timediff((self.start_date,self.start_time),(date,time),24)<0:
            raise KeyError, "Gridded emission file includes (%i,%6.1f) thru (%i,%6.1f); you requested (%i,%6.1f)" % (self.start_date,self.start_time,self.end_date,self.end_time,date,time)
        if chkvar and spc<1 or spc>self.nspec:
            raise KeyError, "Gridded emission file include species 1 thru %i; you requested %i" % (self.nspec,spc)
        
        #self.rffile._newrecord(self.__recordposition(date,time,1,0))
        #start_date,start_time,end_date,end_time=self.rffile.read("ifif")
        self.rffile._newrecord(self.__recordposition(date,time,spc,k))
        
    def read(self):
        """
        Provide direct access to record file read
        """
        return self.rffile.read(self.record_fmt)
        
    def read_into(self,dest):
        """
        Transfer values from current record to dest
        dest - numeric or numpy array
        """
        return read_into(self.rffile,dest,self.id_fmt,self.data_fmt)
        
    def seekandreadinto(self,dest,date=None,time=None,spc=1,k=1):
        """
        see seek and read_into
        """
        self.seek(date,time,spc,k)
        self.read_into(dest)
        
    def seekandread(self,date=None,time=None,spc=1,k=1):
        """
        see seek and read
        """
        self.seek(date,time,spc,k)
        return self.read()

    def itervalues(self):
        for d,t,spc,k in self.__iter__():
            yield self.seekandread(d,t,spc,k)
            
    def iteritems(self):
        for d,t,spc,k in self.__iter__():
            yield d,t,spc,k,self.seekandread(d,t,spc,k)
        
    def iterkeys(self):
        for d,t in self.timerange():
            for spc in range(len(self.spcnames)):
                for k in range(1,self.nlayers+1):
                    yield d,t,spc,k
    __iter__=iterkeys
    def close(self):
        self.rffile.infile.close()
        
    def getArray(self,krange=slice(1,None), nspec=slice(None),nx=slice(None), ny=slice(None)):
        """Method takes slice arguments. Alternatively, takes a hashable object
        with 2 values (e.g., the list: [0,3]). 
        Arguments:
        krange    vertical slice (1 indexed)
        nspec     species  slice (0 indexed)
        nx        column   slice (0 indexed)
        ny        row      slice (0 indexed)
        """
        krange=sliceit(krange)
        nspec=sliceit(nspec)
        nx=sliceit(nx)
        ny=sliceit(ny)
        
        a=zeros(
           (
            self.time_step_count,
            len(xrange(*nspec.indices(self.nspec))),
            len(xrange(*krange.indices(self.nlayers+1))),
            self.ny,
            self.nx)
            ,'f')
        
        for ti,(d,t) in enumerate(self.timerange()):
           for sidx,spc in enumerate(xrange(*nspec.indices(self.nspec))):
               for kidx,k in enumerate(xrange(*krange.indices(self.nlayers+1))):
                    self.seekandreadinto(a[ti,sidx,kidx,...,...],d,t,spc,k)
                  
        return a[...,...,...,ny,nx]

    def timerange(self):
        return timerange((self.start_date,self.start_time),(self.end_date,self.end_time),self.time_step,24)

class uamiv_new(PseudoNetCDFFile):
    def __init__(self,path):
        self.__file=open(path,'r')
        self.__readheader()
        self.__setglobalprops()
        self.__setprivateprops()
        self.__setdimensions()
        self.__settimeprops()
        self.variables=PseudoNetCDFVariables(self.__readonespc,self.__spc_names)
        
    def __readheader(self):
        global_hdr_fmt=dtype([('SPAD1', '>i'), \
                              ('NAME', '>10S4'), \
                              ('NOTE', '>60S4'), \
                              ('IONE', '>i'), \
                              ('NSPEC', '>i'), \
                              ('IBDATE', '>i'), \
                              ('BTIME', '>f'), \
                              ('IEDATE', '>i'), \
                              ('ETIME', '>f'), \
                              ('EPAD1', '>i'), \
                              ('SPAD2', '>i'), \
                              ('RDUM1', '>f'), \
                              ('RDUM2', '>f'), \
                              ('IUTM', '>i'), \
                              ('XORG', '>f'), \
                              ('YORG', '>f'), \
                              ('DELX', '>f'), \
                              ('DELY', '>f'), \
                              ('NX', '>i'), \
                              ('NY', '>i'), \
                              ('NZ', '>i'), \
                              ('IDUM1', '>i'), \
                              ('IDUM2', '>i'), \
                              ('RDUM3', '>f'), \
                              ('RDUM4', '>f'), \
                              ('RDUM5', '>f'), \
                              ('EPAD2', '>i'), \
                              ('SPAD3', '>i'), \
                              ('IONE1', '>i'), \
                              ('IONE2', '>i'), \
                              ('NX2', '>i'), \
                              ('NY2', '>i'), \
                              ('EPAD3', '>i')])
        self.__global_header=fromfile(self.__file,global_hdr_fmt,1)

        spc_fmts=[('SPC%d' % spc, '>10S4') for spc in range(self.__global_header['NSPEC'])]
        spc_hdr_fmt=dtype([('SPAD','>i')]+spc_fmts+[('EPAD','>i')])
        self.__spc_header=fromfile(self.__file,spc_hdr_fmt,1)

    def __setglobalprops(self):
        self.NVARS=nspec=self.__global_header['NSPEC']
        self.NCOLS=nx=self.__global_header['NX']
        self.NROWS=ny=self.__global_header['NY']
        self.NLAYS=nz=max(self.__global_header['NZ'],1)
        self.NAME=self.__global_header['NAME'].reshape(10).view('S1').reshape(10,4)[:,0].tostring()
        self.NOTE=self.__global_header['NOTE'].reshape(60).view('S1').reshape(60,4)[:,0].tostring()
        self.XORIG=self.__global_header['XORG']
        self.YORIG=self.__global_header['YORG']
        self.XCELL=self.__global_header['DELX']
        self.YCELL=self.__global_header['DELY']
        self.SDATE=self.__global_header['IBDATE']+2000000
        self.STIME=self.__global_header['BTIME']*100
        species_name_list=[self.__spc_header['SPC%d' % i].reshape(10).view('S1').reshape(10,4)[:,0].tostring() for i in range(self.NVARS)]
        setattr(self,'VAR-LIST',''.join([spc.ljust(16) for spc in species_name_list]))

    def __setprivateprops(self):
        self.__spc_names=[spc.strip() for spc in array(getattr(self,'VAR-LIST'),ndmin=1).view('S16')]
        self.__spc_ids=dict([(spc,i) for i,spc in enumerate(self.__spc_names)])
        record_buffers=2
        spc_name=10 #10 string4

        self.__date_block=4+record_buffers #4 integer32 + 2 float32
        self.__time_spc_layer_slice=self.NCOLS*self.NROWS+1+spc_name+record_buffers
        self.__time_spc_slice=self.__time_spc_layer_slice*self.NLAYS
        self.__time_slice=self.__time_spc_slice*self.NVARS
        
        self.__time_block=self.__date_block+self.__time_slice
        
        self.__inc=(self.__time_block-self.__time_spc_slice)*4
        self.__data_start=self.__file.tell()
        self.__file.seek(0,2)
        self.__flen=self.__file.tell()

    def __settimeprops(self):
        ntimes=(self.__flen-self.__data_start)/4./self.__time_block
        self.createDimension('TSTEP',ntimes)
        self.NSTEPS=ntimes
    
    def __spcstart(self,nspc):
        return self.__data_start+(self.__date_block+nspc*self.__time_spc_slice)*4

    def __seektospc(self,spc):
        start=self.__spcstart(self.__spc_ids[spc.strip()])
        self.__file.seek(start,0)
    
    def __readonespc(self,spc):
        self.__file.seek(0,2)
        print self.__file.tell()
        self.__seektospc(spc)
        print self.__file.tell()
        vals=zeros((self.NSTEPS,self.__time_spc_slice),'>f')
        for hri in range(self.NSTEPS):
            vals[hri,:]=fromfile(self.__file,'>f',self.__time_spc_slice)
            self.__file.seek(self.__inc,1)
        
        vals=vals.reshape(self.NSTEPS,self.NLAYS,self.__time_spc_layer_slice)
        vals=vals[:,:,12:-1]
        vals=vals.reshape(self.NSTEPS,self.NLAYS,self.NROWS,self.NCOLS)
        units={'AVERAGE   ':'ppm','AIRQUALITY':'ppm','EMISSIONS ':'mol'}[self.NAME].ljust(16)
        vals=PseudoIOAPIVariable(self,spc,'f',('TSTEP','LAY','ROW','COL'),values=vals, units = units)

        return vals

    def __setdimensions(self):
        self.dimensions={}
        self.createDimension('LAY',self.NLAYS)
        self.createDimension('ROW',self.NROWS)
        self.createDimension('COL',self.NCOLS)

    

        

class TestuamivRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testGE(self):
        emissfile=uamiv('../../../../testdata/ei/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
        self.assert_((emissfile.variables['NO'].mean(1).mean(1).mean(1)==array([  52.05988312,   51.58646774,   51.28796387,   55.63090134,
         63.95315933,  105.3456192 ,  158.26776123,  152.04057312,
         147.32403564,  154.80661011,  164.03274536,  171.88658142,
         174.36567688,  180.03359985,  173.81938171,  180.50257874,
         178.56637573,  161.35736084,  110.38669586,   97.90225983,
         89.08138275,   81.10474396,   73.36611938,   58.82622528],dtype='f')).all())
    def testGENew(self):
        emissfile=uamiv_new('../../../../testdata/ei/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
        self.assert_((emissfile.variables['NO'].mean(1).mean(1).mean(1)==array([  52.05988312,   51.58646774,   51.28796387,   55.63090134,
         63.95315933,  105.3456192 ,  158.26776123,  152.04057312,
         147.32403564,  154.80661011,  164.03274536,  171.88658142,
         174.36567688,  180.03359985,  173.81938171,  180.50257874,
         178.56637573,  161.35736084,  110.38669586,   97.90225983,
         89.08138275,   81.10474396,   73.36611938,   58.82622528],dtype='f')).all())
       
if __name__ == '__main__':
    unittest.main()
