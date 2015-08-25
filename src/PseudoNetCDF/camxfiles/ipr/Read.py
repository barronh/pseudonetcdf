__all__=['ipr']
__doc__ = """
.. _Read
:mod:`Read` -- ipr Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              ipr files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from numpy import zeros,array,where,memmap,newaxis,dtype,fromfile,dtype

#This Package modules
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.util import cartesian
from PseudoNetCDF.camxfiles.units import get_uamiv_units
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class ipr(PseudoNetCDFFile):
    """
    ipr provides a PseudoNetCDF interface for CAMx
    ipr files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> ipr_path = 'camx_ipr.bin'
        >>> rows,cols = 65,83
        >>> iprfile = ipr(ipr_path,rows,cols)
        >>> iprfile.variables.keys()
        ['TFLAG', 'SPAD_O3', 'DATE_O3', 'TIME_O3', 'SPC_O3', 
         'PAGRID_O3', 'NEST_O3', 'I_O3', 'J_O3', 'K_O3', 
         'INIT_O3', 'CHEM_O3', 'EMIS_O3', 'PTEMIS_O3', 
         'PIG_O3', 'WADV_O3', 'EADV_O3', 'SADV_O3', 'NADV_O3', 
         'BADV_O3', 'TADV_O3', 'DIL_O3', 'WDIF_O3', 'EDIF_O3', 
         'SDIF_O3', 'NDIF_O3', 'BDIF_O3', 'TDIF_O3', 'DDEP_O3', 
         'WDEP_O3', 'INORGACHEM_O3', 'ORGACHEM_O3', 'AQACHEM_O3', 
         'FCONC_O3', 'UCNV_O3', 'AVOL_O3', 'EPAD_O3']
        >>> v = iprfile.variables['CHEM_O3']
        >>> tflag = iprfile.variables['TFLAG']
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
        >>> iprfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    __ipr_record_type ={
        24: dtype(
              dict(
                 names=['SPAD', 'DATE', 'TIME', 'SPC', 'PAGRID', 'NEST', 'I', 'J', 'K', 'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'WADV', 'EADV', 'SADV', 'NADV', 'BADV', 'TADV', 'DIL', 'WDIF', 'EDIF', 'SDIF', 'NDIF', 'BDIF', 'TDIF', 'DDEP', 'WDEP', 'AERCHEM', 'FCONC', 'UCNV', 'AVOL', 'EPAD'], 
                 formats=['>i', '>i', '>f', '>S10', '>i', '>i', '>i', '>i', '>i', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>i']
                )
              ), 
        26: dtype(
              dict(
                 names=['SPAD', 'DATE', 'TIME', 'SPC', 'PAGRID', 'NEST', 'I', 'J',  'K', 'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'WADV', 'EADV', 'SADV', 'NADV', 'BADV', 'TADV', 'DIL', 'WDIF', 'EDIF', 'SDIF', 'NDIF', 'BDIF', 'TDIF', 'DDEP', 'WDEP', 'INORGACHEM', 'ORGACHEM', 'AQACHEM', 'FCONC', 'UCNV', 'AVOL', 'EPAD'], 
                 formats=['>i', '>i', '>f', '>S10', '>i', '>i', '>i', '>i', '>i', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>i']
                )
              )
         }

    def __init__(self,rf,proc_dict=None,units='umol/m**3', oldnames = False, **props):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        
        Keywords (i.e., props) for projection: P_ALP, P_BET, P_GAM, XCENT, YCENT, XORIG, YORIG, XCELL, YCELL
        """
        if proc_dict!=None:
            self.proc_dict=proc_dict
        else:
            self.proc_dict=None
        self.__rffile=open(rf, 'rb')
        self.__rffile.seek(0,2)
        if self.__rffile.tell()<2147483648L:
            warn("For greater speed on files <2GB use ipr_memmap")
        self.__rffile.seek(0,0)
        self.units=units
        self.__readheader()
        self.__setDomain__()
        self.__gettimestep()
        tdim = self.createDimension('TSTEP',self.NSTEPS)
        tdim.setunlimited(True)
        self.createDimension('DATE-TIME',2)
        self.createDimension('VAR',self.NSPCS*self.NPROCESS)
        varkeys=["_".join([j[1],j[0]]) for j in [i for i in cartesian([i.strip() for i in self.spcnames['SPECIES'].tolist()],self.proc_dict.keys())]]+['TFLAG']
        self.variables=PseudoNetCDFVariables(self.__variables,varkeys)
        for k, v in props.iteritems():
            setattr(self, k, v)
        add_cf_from_ioapi(self)
        

    def __variables(self,proc_spc):
        if proc_spc=='TFLAG':
            time=self.variables['TIME_%s'  % self.spcnames[0][1].strip()]
            date=self.variables['DATE_%s'  % self.spcnames[0][1].strip()]
            self.variables['TFLAG']=PseudoNetCDFVariable(self,'proc_spc','i',('TSTEP','VAR','DATE-TIME'),values=ConvertCAMxTime(date[:,0,0,0],time[:,0,0,0],len(self.dimensions['VAR'])))
            return self.variables['TFLAG']
        
        self.variables.clear()

        for k in self.proc_dict:
            proc=proc_spc[:len(k)]
            spc=proc_spc[len(k)+1:]
            if proc==k and spc.ljust(10) in self.spcnames['SPECIES'].tolist():
                spcprocs=self.__readalltime(spc)
                for p,plong in self.proc_dict.iteritems():
                    var_name=p+'_'+spc
                    # IPR units are consistent with 'IPR'
                    if p == 'UCNV':
                        units = 'm**3/mol'
                    elif p == 'AVOL':
                        units = 'm**3'
                    else:
                        units = get_uamiv_units('IPR', spc)
                    self.variables[var_name] = PseudoNetCDFVariable(self,var_name,'f',('TSTEP','LAY','ROW','COL'),
                                              values=spcprocs[p],
                                              units=units,
                                              var_desc=(var_name).ljust(16),
                                              long_name=(var_name).ljust(16))
                del spcprocs
                return self.variables[proc_spc]
        raise KeyError, "Bad!"

    def __readonetime(self,ntime,spc):
        newpos = self.__start(ntime,spc)
        oldpos = self.__rffile.tell()
        relmove = newpos - oldpos
        self.__rffile.seek(relmove, 1)
        return fromfile(self.__rffile,dtype=self.__ipr_record_type,count=self.__block3d).reshape(self.NROWS,self.NCOLS,self.NLAYS).swapaxes(0,2).swapaxes(1,2)
    
    def __readalltime(self,spc):
        out=zeros((self.NSTEPS, self.NLAYS, self.NROWS, self.NCOLS), dtype = self.__ipr_record_type)
        for it in range(self.NSTEPS):
            out[it] = self.__readonetime(it,spc)
        return out

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
                'WADV': 'West boundary advection', 
                'EADV': 'East boundary advection', 
                'SADV': 'South boundary advection', 
                'NADV': 'North boundary advection', 
                'BADV': 'Bottom boundary advection', 
                'TADV': 'Top boundary advection', 
                'DIL': 'Dilution in the vertical', 
                'WDIF': 'West boundary diffusion', 
                'EDIF': 'East boundary diffusion', 
                'SDIF': 'South boundary diffusion', 
                'NDIF': 'North boundary diffusion', 
                'BDIF': 'Bottom boundary diffusion', 
                'TDIF': 'Top boundary diffusion', 
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
            elif self.NPROCESS[0] == 26:
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
        self.NCOLS=len(self.dimensions['COL'])
        self.NROWS=len(self.dimensions['ROW'])
        self.NLAYS=len(self.dimensions['LAY'])
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
        temp=fromfile(self.__rffile,dtype=self.__ipr_record_type,count=len(self.dimensions['LAY'])*len(self.dimensions['ROW'])*len(self.dimensions['COL'])+1)
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
