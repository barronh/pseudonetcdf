__all__=['ipr']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- ipr Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              ipr files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
import unittest
import struct
from warnings import warn
from mmap import error as MemmapLimitError

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.util import cartesian
from PseudoNetCDF.camxfiles.units import get_uamiv_units
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
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
        >>> iprfile = ipr(ipr_path)
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
    
    id_fmt="if10s5i"
    dt_fmt="if"
    data_fmt="f"
    
    def __init__(self,rf,multi=False, **props):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info

        Keywords (i.e., props) for projection: P_ALP, P_BET, P_GAM, XCENT, YCENT, XORIG, YORIG, XCELL, YCELL
        """
        self.__rffile=OpenRecordFile(rf)
        self.__readheader()
        self.__ipr_record_type={
            24: dtype(
                        dict(
                            names=['SPAD', 'DATE', 'TIME', 'SPC', 'PAGRID', 'NEST', 'I', 'J', 'K', 
                                    'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'WADV', 'EADV', 'SADV', 
                                    'NADV', 'BADV', 'TADV', 'DIL', 'WDIF', 'EDIF', 'SDIF', 'NDIF', 
                                    'BDIF', 'TDIF', 'DDEP', 'WDEP', 'AERCHEM', 'FCONC', 'UCNV', 'AVOL', 
                                    'EPAD'], 
                            formats=['>i', '>i', '>f', '>S10', '>i', '>i', '>i', '>i', '>i', 
                                    '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f',
                                    '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f',
                                    '>f', '>f', '>f', '>f', '>f', '>f', '>i'])),
            26: dtype(
                        dict(
                            names=['SPAD', 'DATE', 'TIME', 'SPC', 'PAGRID', 'NEST', 'I', 'J', 'K', 
                                    'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'WADV', 'EADV', 'SADV', 
                                    'NADV', 'BADV', 'TADV', 'DIL', 'WDIF', 'EDIF', 'SDIF', 'NDIF', 
                                    'BDIF', 'TDIF', 'DDEP', 'WDEP', 'INORGACHEM', 'ORGACHEM', 'AQACHEM', 'FCONC', 'UCNV', 'AVOL', 
                                    'EPAD'], 
                            formats=['>i', '>i', '>f', '>S10', '>i', '>i', '>i', '>i', '>i', 
                                    '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f',
                                    '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f',
                                    '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>f', '>i']))
                                }[len(self.prcnames)]

        prcs=['SPAD', 'DATE', 'TIME', 'PAGRID', 'NEST', 'I', 'J', 'K', 
                'INIT', 'CHEM', 'EMIS', 'PTEMIS', 'PIG', 'WADV', 'EADV', 'SADV', 
                'NADV', 'BADV', 'TADV', 'DIL', 'WDIF', 'EDIF', 'SDIF', 'NDIF', 
                'BDIF', 'TDIF', 'DDEP', 'WDEP']+{24: ['AERCHEM'], 26: ['INORGACHEM', 'ORGACHEM', 'AQACHEM']}[len(self.prcnames)]+['FCONC', 'UCNV', 'AVOL', 
                'EPAD']
        varkeys=['_'.join(i) for i in cartesian(prcs,self.spcnames)]
        varkeys+=['SPAD','DATE','TIME','PAGRID','NEST','I','J','K','TFLAG']
        self.groups = {}
        NSTEPS = len([i_ for i_ in self.timerange()])
        NVARS = len(varkeys)
        self.createDimension('VAR', NVARS)
        self.createDimension('DATE-TIME', 2)
        self.createDimension('TSTEP', NSTEPS)
        padatatype = []
        pavarkeys = []
        for di, domain in enumerate(self.padomains):
            dk = 'PA%02d' % di
            prefix = dk + '_'
            grp = self.groups[dk] = PseudoNetCDFFile()
            pavarkeys.extend([prefix + k for k in varkeys])
            grp.createDimension('VAR', NVARS)
            grp.createDimension('DATE-TIME', 2)
            grp.createDimension('TSTEP', NSTEPS)
            grp.createDimension('COL', domain['iend'] - domain['istart'] + 1)
            grp.createDimension('ROW', domain['jend'] - domain['jstart'] + 1)
            grp.createDimension('LAY', domain['tlay'] - domain['blay'] + 1)
            padatatype.append((dk, self.__ipr_record_type, (len(grp.dimensions['ROW']), len(grp.dimensions['COL']), len(grp.dimensions['LAY']))))
            if len(self.padomains) == 1:
                self.createDimension('COL', domain['iend']-domain['istart']+1)
                self.createDimension('ROW', domain['jend']-domain['jstart']+1)
                self.createDimension('LAY', domain['tlay']-domain['blay']+1)
            exec("""def varget(k):
                return self._ipr__variables('%s', k)""" % dk, dict(self = self), locals())
            if len(self.padomains) == 1:
                self.variables = PseudoNetCDFVariables(varget,varkeys)
            else:
                grp.variables = PseudoNetCDFVariables(varget,varkeys)
        
        self.__memmaps=memmap(self.__rffile.infile.name,dtype(padatatype),'r',self.data_start_byte).reshape(NSTEPS, len(self.spcnames))
        for k, v in props.iteritems():
            setattr(self, k, v)
        try:
            add_cf_from_ioapi(self)
        except:
            pass

    def __del__(self):
        try:
            self.__memmaps.close()
            del self.__memmaps
        except:
            pass

    def __decorator(self,name,pncfv):
        spc = name.split('_')[-1]
        prc = name.split('_')[0]
        # IPR units are consistent with 'IPR'
        if prc == 'UCNV':
            units = 'm**3/mol'
        elif prc == 'AVOL':
            units = 'm**3'
        else:
            units = get_uamiv_units('IPR', spc)
        decor=lambda k: dict(units=units, var_desc=k.ljust(16), long_name=k.ljust(16))
        for k,v in decor(name).iteritems():
            setattr(pncfv,k,v)        
        return pncfv
        
    def __variables(self,pk, proc_spc):
        if proc_spc in self.__ipr_record_type.names:
            proc=proc_spc
            proc_spc=proc_spc+'_'+self.spcnames[0]
            return PseudoNetCDFVariable(self,proc_spc,'f',('TSTEP','LAY','ROW','COL'),values=self.__memmaps[pk][:,0,:,:,:][proc].swapaxes(1, 3).swapaxes(2, 3))
        if proc_spc=='TFLAG':
            thisdate = self.__memmaps[pk][:,0,:,:,:]['DATE'].swapaxes(1, 3).swapaxes(2, 3)[..., 0, 0, 0]
            thistime = self.__memmaps[pk][:,0,:,:,:]['TIME'].swapaxes(1, 3).swapaxes(2, 3)[..., 0, 0, 0]
            return ConvertCAMxTime(thisdate, thistime, len(self.groups[pk].dimensions['VAR']))
        for k in self.__ipr_record_type.names:
            proc=proc_spc[:len(k)]
            spc=proc_spc[len(k)+1:]
            if proc==k and spc in self.spcnames:
                spc=self.spcnames.index(spc)
                dvals = self.__memmaps[pk][:,spc][proc].swapaxes(1, 3).swapaxes(2, 3)
                return self.__decorator(proc_spc,PseudoNetCDFVariable(self,proc_spc,'f',('TSTEP','LAY','ROW','COL'),values=dvals))
        raise KeyError, "Bad!"
                
                
    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        
        self.runmessage=self.__rffile.read("80s")
        self.start_date,self.start_time,self.end_date,self.end_time=self.__rffile.read("ifif")
        
        self.grids=[]
        for grid in range(self.__rffile.read("i")[-1]):
            self.grids.append(
                            dict(
                                zip(
                                    ['orgx','orgy','ncol','nrow','xsize','ysize'], 
                                    self.__rffile.read("iiiiii")
                                    )
                                )
                            )
        
        self.spcnames = []
        for spc in range(self.__rffile.read("i")[-1]):
            self.spcnames.append(self.__rffile.read("10s")[-1].strip())
            
        self.nspec=len(self.spcnames)
        self.padomains=[]
        
        for padomain in range(self.__rffile.read("i")[-1]):
            self.padomains.append(
                                dict(
                                    zip(
                                        ['grid','istart','iend','jstart','jend','blay','tlay'],
                                        self.__rffile.read("iiiiiii")
                                        )
                                    )
                                )
        self.activedomain=self.padomains[0]
        self.prcnames=[]
        
        for i in range(self.__rffile.read('i')[-1]):
            self.prcnames.append(self.__rffile.read('25s')[-1].strip())
        
        self.data_start_byte=self.__rffile.record_start
        self.record_fmt=self.id_fmt + str(len(self.prcnames)) + self.data_fmt
        self.record_size=self.__rffile.record_size
        self.SDATE,self.STIME,dummy,dummy,dummy,dummy,dummy,dummy=self.__rffile.read(self.id_fmt)
        self.__rffile.previous()
        self.TSTEP=100.
        self.padded_size=self.record_size+8
        domain=self.padomains[0]
        self.records_per_time=self.nspec*(domain['iend']-domain['istart']+1)*(domain['jend']-domain['jstart']+1)*(domain['tlay']-domain['blay']+1)
        self.time_data_block=self.padded_size*self.records_per_time
        self.time_step=100.

    def timerange(self):
        return timerange((self.start_date,self.start_time+self.time_step),timeadd((self.end_date,self.end_time),(0,self.time_step)),self.time_step)

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testIPR(self):
        warn('Test not implemented')
       
if __name__ == '__main__':
    unittest.main()
