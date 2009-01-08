__all__=['ipr']

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
from numpy import zeros,array,where,memmap,newaxis,dtype,fromfile,dtype
try:
    from Scientific.IO.NetCDF import NetCDFFile as ncf
except:
    from pynetcdf import NetCDFFile as ncf

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd,timerange
from pyPA.utils.util import cartesian,sliceit
from pyPA.utils.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class ipr(PseudoNetCDFFile):
    __ipr_record_type={
        24: dtype(
              dict(
                 names=['SPAD', 'DATE', 'TIME', 'SPC', 'PAGRID', 'NEST', 'I', 'J', 'K', 'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'A_W', 'A_E', 'A_S', 'A_N', 'A_B', 'A_T', 'DIL', 'D_W', 'D_E', 'D_S', 'D_N', 'D_B', 'D_T', 'DDEP', 'WDEP', 'AERCHEM', 'FCONC', 'UCNV', 'AVOL', 'EPAD'], 
                 formats=['>i', '>i', '>f', '>S10', '>i', '>i', '>i', '>i', '>i', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>i']
                )
              ), 
        26: dtype(
              dict(
                 names=['SPAD', 'DATE', 'TIME', 'SPC', 'PAGRID', 'NEST', 'I', 'J',  'K', 'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'A_W', 'A_E', 'A_S', 'A_N', 'A_B', 'A_T', 'DIL', 'D_W', 'D_E', 'D_S', 'D_N', 'D_B', 'D_T', 'DDEP', 'WDEP', 'INORGACHEM', 'ORGACHEM', 'AQACHEM', 'FCONC', 'UCNV', 'AVOL', 'EPAD'], 
                 formats=['>i', '>i', '>f', '>S10', '>i', '>i', '>i', '>i', '>i', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>i']
                )
              )
         }
    def __init__(self,rf,proc_dict=None,units='umol/hr'):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        if proc_dict!=None:
            self.proc_dict=proc_dict
        else:
            self.proc_dict=None
        self.__rffile=file(rf)
        self.__rffile.seek(0,2)
        if self.__rffile.tell()<2147483648L:
            warn("For greater speed on files <2GB use ipr_memmap")
        self.__rffile.seek(0,0)
        self.units=units
        self.dimensions={}
        self.__readheader()
        self.__setDomain__()
        self.__gettimestep()
        self.createDimension('TSTEP',self.NSTEPS)
        self.createDimension('DATE-TIME',2)
        self.createDimension('VAR',self.NSPCS*self.NPROCESS)
        varkeys=["_".join([j[1],j[0]]) for j in [i for i in cartesian([i.strip() for i in self.spcnames['SPECIES'].tolist()],self.proc_dict.keys())]]+['TFLAG']
        self.__invariables=[]
        self.variables=PseudoNetCDFVariables(self.__variables,varkeys)
        

    def __variables(self,proc_spc):
        if proc_spc=='TFLAG':
            time=self.variables['TIME_%s'  % self.spcnames[0][1].strip()]
            date=self.variables['DATE_%s'  % self.spcnames[0][1].strip()]
            self.variables['TFLAG']=PseudoNetCDFVariable(self,'proc_spc','i',('TSTEP','VAR','DATE-TIME'),values=ConvertCAMxTime(date[:,0,0,0],time[:,0,0,0],self.dimensions['VAR']))
            return self.variables['TFLAG']
        for k in self.__invariables:
            try:
                del self.variables[k]
            except:
                pass
                
        self.__invariables=[]
        for k in self.proc_dict:
            proc=proc_spc[:len(k)]
            spc=proc_spc[len(k)+1:]
            if proc==k and spc.ljust(10) in self.spcnames['SPECIES'].tolist():
                spcprocs=self.__readalltime(spc)
                for p,plong in self.proc_dict.iteritems():
                    var_name=p+'_'+spc
                    tmpv=PseudoNetCDFVariable(self,var_name,'f',('TSTEP','LAY','ROW','COL'),values=spcprocs[p])
                    self.__invariables.append(var_name)
                    tmpv.units='umol/hr'
                    tmpv.var_desc=(var_name).ljust(16)
                    tmpv.long_name=(var_name).ljust(16)
                    self.variables[var_name]=tmpv
                    tmpv=None
                del spcprocs
                return self.variables[proc_spc]
        raise KeyError, "Bad!"

    def __readonetime(self,ntime,spc):
        self.__rffile.seek(self.__start(ntime,spc),0)
        return fromfile(self.__rffile,dtype=self.__ipr_record_type,count=self.__block3d).reshape(self.NROWS,self.NCOLS,self.NLAYS).swapaxes(0,2).swapaxes(1,2)
    
    def __readalltime(self,spc):
        out=[]
        for it in range(self.NSTEPS):
            out.append(self.__readonetime(it,spc))
        return array(out,dtype=self.__ipr_record_type)

    def __start(self,ntime,spc):
        nspec=self.spcnames['SPECIES'].tolist().index(spc.ljust(10))
        return self.__data_start_byte+(long(ntime)*self.__block4d+self.__block3d*nspec)*self.__ipr_record_type.itemsize
    
    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        
        self.runmessage=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','RUNMESSAGE','EPAD'],formats=['>i','>80S','>i'])),count=1)['RUNMESSAGE']
        dates=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','SDATE','STIME','EDATE','ETIME','EPAD'],formats=['>i','>i','>f','>i','>f','>i'])),count=1)
        self.SDATE=dates['SDATE']+2000000
        self.STIME=dates['STIME']
        self.EDATE=dates['EDATE']+2000000
        self.ETIME=dates['ETIME']
        
        self.__grids=[]
        self.NGRIDS=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','NGRIDS','EPAD'],formats=['>i']*3)),count=1)['NGRIDS']
        for grid in range(self.NGRIDS):
            self.__grids.append( 
                                    fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','orgx','orgy','ncol','nrow','xsize','ysize','EPAD'],formats=['>i','>i','>i','>i','>i','>i','>i','>i'])),count=1)
                            )
        
        self.spcnames = []
        self.NSPCS=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','NSPCS','EPAD'],formats=['>i','>i','>i'])),count=1)['NSPCS'][0]
        self.spcnames=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','SPECIES','EPAD'],formats=['>i','>10S','>i'])),count=self.NSPCS)
        
        self.padomains=[]
        self.NPADOMAINS=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','NPADOMAINS','EPAD'],formats=['>i','>i','>i'])),count=1)['NPADOMAINS'][0]
        self.__padomains=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','grid','istart','iend','jstart','jend','blay','tlay','EPAD'],formats=['>i','>i','>i','>i','>i','>i','>i','>i','>i'])),count=self.NPADOMAINS)
        self.__activedomain=self.__padomains[0]
        self.prcnames=[]
        self.NPROCESS=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','NPRCS','EPAD'],formats=['>i','>i','>i'])),count=1)['NPRCS']
        self.__ipr_record_type=self.__ipr_record_type[self.NPROCESS[0]]
        
        if self.proc_dict is None:
            self.proc_dict={
                'INIT': 'Initial concentration', 
                'CHEM': 'Chemistry', 
                'EMIS': 'Area emissions',
                'PTEMIS': 'Point source emissions', 
                'PIG': 'Plume-in-Grid change', 
                'A_W': 'West boundary advection', 
                'A_E': 'East boundary advection', 
                'A_S': 'South boundary advection', 
                'A_N': 'North boundary advection', 
                'A_B': 'Bottom boundary advection', 
                'A_T': 'Top boundary advection', 
                'DIL': 'Dilution in the vertical', 
                'D_W': 'West boundary diffusion', 
                'D_E': 'East boundary diffusion', 
                'D_S': 'South boundary diffusion', 
                'D_N': 'North boundary diffusion', 
                'D_B': 'Bottom boundary diffusion', 
                'D_T': 'Top boundary diffusion', 
                'DDEP': 'Dry deposition', 
                'WDEP': 'Wet deposition', 
                'FCONC': 'Final concentration', 
                'UCNV': 'Units conversion', 
                'AVOL': 'Average cell volume', 
                'DATE': 'DATE', 
                'TIME': 'TIME', 
                'K': 'K', 
                'J': 'J', 
                'I': 'I'
                }
            if self.NPROCESS[0] == 24:
                self.proc_dict['AERCHEM'] = 'Aerosol chemistry'
            elif self.NPROCESS[0] == 24:
                self.proc_dict['INORGACHEM'] = 'Inorganic Aerosol chemistry'
                self.proc_dict['ORGACHEM'] = 'Organic Aerosol chemistry'
                self.proc_dict['AQACHEM'] = 'Aqueous Aerosol chemistry'
            else:
                warn('Unknown version; cannot add aerosol chemistry') 

        
        self.prcnames=fromfile(self.__rffile,dtype=dtype(dict(names=['SPAD','PROCESS','EPAD'],formats=['>i','>25S','>i'])),count=self.NPROCESS)
        self.__data_start_byte=self.__rffile.tell()
        
    def __setDomain__(self,id=0):
        self.__activedomain=self.__padomains[id]
        self.createDimension('COL',self.__activedomain['iend']-self.__activedomain['istart']+1)
        self.createDimension('ROW',self.__activedomain['jend']-self.__activedomain['jstart']+1)
        self.createDimension('LAY',self.__activedomain['tlay']-self.__activedomain['blay']+1)
        self.NCOLS=self.dimensions['COL']
        self.NROWS=self.dimensions['ROW']
        self.NLAYS=self.dimensions['LAY']
        self.__block3d=self.NLAYS*self.NROWS*self.NCOLS
        self.__block4d=self.__block3d*self.NSPCS
        
    def __gettimestep(self):
        """
        Header information provides start and end date, but does not
        indicate the increment between.  This routine reads the first
        and second date/time and initializes variables indicating the
        timestep length and the anticipated number.
        """
        self.__rffile.seek(self.__data_start_byte,0)
        temp=fromfile(self.__rffile,dtype=self.__ipr_record_type,count=self.dimensions['LAY']*self.dimensions['ROW']*self.dimensions['COL']+1)
        self.TSTEP=timediff((self.SDATE,self.STIME),(temp[-1]['DATE']+2000000,temp[-1]['TIME']))
        self.NSTEPS=int(timediff((self.SDATE,self.STIME),(self.EDATE,self.ETIME))/self.TSTEP)

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testIPR(self):
        emissfile=ipr('../../../../testdata/ei/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
        self.assert_(1==2)
        
if __name__ == '__main__':
    unittest.main()
