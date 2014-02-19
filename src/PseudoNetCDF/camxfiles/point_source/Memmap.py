__all__=['point_source']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- point_source Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              point_source files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan
import numpy as np
#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

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
    
    __emiss_hdr_fmt=dtype(dict(names=['SPAD','name','note','itzon','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>(10,4)S1','>(60,4)S1','>i','>i','>i','>f','>i','>f','>i']))
    __grid_hdr_fmt=dtype(dict(names=['SPAD','plon','plat','iutm','xorg','yorg','delx','dely','nx','ny','nz', 'iproj','istag','tlat1','tlat2','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))
    __cell_hdr_fmt=dtype(dict(names=['SPAD','ione1','ione2','nx','ny','EPAD'],formats=['>i','>i','>i','>i','>i','>i']))
    
    __spc_hdr_fmt=dtype('(10,4)>S1')
    __time_hdr_fmt=dtype(dict(names=['SPAD','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))
    __nstk_hdr_fmt=dtype(dict(names=['SPAD','ione','nstk','EPAD'],formats=['>i','>i','>i','>i']))
    
    
    __stk_prop_fmt=dtype(dict(names=['XSTK','YSTK','HSTK','DSTK','TSTK','VSTK'],formats=['>f','>f','>f','>f','>f','>f']))
    
    __stk_time_prop_fmt=dtype(dict(names=['IONE','ITWO','KCELL','FLOW','PLMHT'],formats=['>i','>i','>f','>f','>f']))

    def __init__(self,rf):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.__rffile=rf
        
        
        self.variables={}
        self.__memmap=memmap(self.__rffile,'>f','r')
        self.__globalheader()
        varkeys=['ETFLAG', 'TFLAG', 'XSTK', 'YSTK', 'HSTK', 'DSTK', 'TSTK', 'VSTK', 'IONE', 'ITWO', 'KCELL', 'FLOW', 'PLMHT', 'NSTKS']+[i.strip() for i in self.__spc_names]
        self.SPC_NAMES = ''.join([s.ljust(16) for s in self.__spc_names])
        setattr(self,'VAR-LIST',"".join([i.ljust(16) for i in varkeys]))
        self.variables=PseudoNetCDFVariables(self.__variables,varkeys)
        self.__time_stks()
        self.createDimension('NSTK',self.__nstk_hdr['nstk'])
        
        
    
    def __globalheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        offset=0

        self.__emiss_hdr=self.__memmap[offset:offset+self.__emiss_hdr_fmt.itemsize/4].view(self.__emiss_hdr_fmt)
        offset+=self.__emiss_hdr.nbytes/4

        self.__grid_hdr=self.__memmap[offset:offset+self.__grid_hdr_fmt.itemsize/4].view(self.__grid_hdr_fmt)
        offset+=self.__grid_hdr.nbytes/4
        self.NAME = self.__emiss_hdr['name'][0, :, 0].copy().view('S10')[0]
        self.NOTE = self.__emiss_hdr['note'][0, :, 0].copy().view('S60')[0]
        self.XORIG=self.__grid_hdr['xorg'][0]
        self.YORIG=self.__grid_hdr['yorg'][0]
        self.XCELL=self.__grid_hdr['delx'][0]
        self.YCELL=self.__grid_hdr['dely'][0]
        self.NCOLS = self.__grid_hdr['nx'][0]
        self.NROWS = self.__grid_hdr['ny'][0]
        self.NLAYS = self.__grid_hdr['nz'][0]
        self.PLON = plon = self.__grid_hdr['plon'][0]
        self.PLAT = plat = self.__grid_hdr['plat'][0]
        self.TLAT1 = tlat1 = self.__grid_hdr['tlat1'][0]
        self.TLAT2 = tlat2 = self.__grid_hdr['tlat2'][0]
        self.IUTM = iutm = self.__grid_hdr['iutm'][0]
        self.ISTAG = istag = self.__grid_hdr['istag'][0]        
        self.CPROJ = cproj = self.__grid_hdr['iproj'][0]
        self.__cell_hdr=self.__memmap[offset:offset+self.__cell_hdr_fmt.itemsize/4].view(self.__cell_hdr_fmt)
        offset+=self.__cell_hdr.nbytes/4+1

        nspec=self.__emiss_hdr['nspec'][0]
        self.ITZON=self.__emiss_hdr['itzon'][0]


        self.__spc_hdr=self.__memmap[offset:offset+self.__spc_hdr_fmt.itemsize/4*nspec].view('>S1').reshape(nspec, 10, 4)
        offset+=self.__spc_hdr.nbytes/4+1

        self.__spc_names=[str(np.char.strip(spc[:,0].copy().view('S10')[0])) for spc in self.__spc_hdr]
        self.__nstk_hdr=self.__memmap[offset:offset+self.__nstk_hdr_fmt.itemsize/4].view(self.__nstk_hdr_fmt)
        offset+=self.__nstk_hdr.nbytes/4+1

        nstk=self.__nstk_hdr['nstk']
        self.createDimension('NSTK', nstk)

        self.__nstk_hdr['nstk']
        self.__stk_props=self.__memmap[offset:offset+self.__stk_prop_fmt.itemsize/4*nstk].reshape(nstk,self.__stk_prop_fmt.itemsize/4).view(self.__stk_prop_fmt)
        offset+=self.__stk_props.nbytes/4+1

        self.__data_start=offset
        
        self.SDATE=self.__emiss_hdr['ibdate']
        self.STIME=self.__emiss_hdr['btime']
        self.EDATE=self.__emiss_hdr['iedate']
        self.ETIME=self.__emiss_hdr['etime']
        
    def __getspcidx(self,spc):
        return self.__spc_names.index(spc)

    def __time_stks(self):
        i=offset=0
        nspcs=len(self.__spc_names)
        nstks=len(self.dimensions['NSTK'])
        date_block_size=6
        stk_block_size=4
        stk_props_size=2+nstks*5
        emiss_block_size=nspcs*(nstks+13)
        hour_block_size=date_block_size+stk_block_size+stk_props_size+emiss_block_size
        data=self.__memmap[self.__data_start:]
        data=data.reshape(data.size/hour_block_size,hour_block_size)
        ntimes=data.shape[0]
        self.createDimension('TSTEP',ntimes)
        self.createDimension('DATE-TIME', 2)
        start=0
        end=date_block_size
        date_times=data[:,start:end]
        dates=date_times[:,[1,3]].view('>i')
        times=date_times[:,[2,4]]
        
        start=end
        end=start+stk_block_size
        nstk_hdr=data[:,start:end]
        if not (nstks==nstk_hdr[:,2:3].view('>i')).all():
            raise ValueError, "Number of stacks varies with time"
        start=end
        end=start+stk_props_size
        self.__hourly_stk_props=data[:,start:end][:,1:-1].reshape(ntimes,nstks,5)
        
        start=end
        end=start+emiss_block_size
        if not end==data.shape[1]:
            raise ValueError, "Incorrect shape"
        self.__emiss_data=data[:,start:].reshape(ntimes,nspcs,13+nstks)[:,:,12:-1]
        bdates=dates[:,0]
        btimes=times[:,0]
        edates=dates[:,1]
        etimes=times[:,1]
        self.NSTEPS=ntimes
        
        self.createDimension('VAR',len(self.__spc_names)+3)
        
        self.variables['TFLAG']=ConvertCAMxTime(bdates,btimes,len(self.dimensions['VAR']))
        self.variables['ETFLAG']=ConvertCAMxTime(edates,etimes,len(self.dimensions['VAR']))
        v=self.variables['NSTKS']=PseudoNetCDFVariable(self,'NSTKS','i',('TSTEP',),values=array(nstk_hdr))
        v.units='#'.ljust(16)
        v.long_name='NSTKS'.ljust(16)
        v.var_desc=v.long_name
    
    def __variables(self,k):
        if k in ['TFLAG','ETFLAG','NSTKS']:
            return self.variables[k]
        elif k in ['XSTK','YSTK','HSTK','DSTK','TSTK','VSTK']:
            v=PseudoNetCDFVariable(self,k,'f',('NSTK',),values=self.__stk_props[k].ravel())
            v.units={'XSTK':'m','YSTK':'m','HSTK':'m','DSTK':'m','TSTK':'K','VSTK':'m/h'}[k]
            v.long_name=k.ljust(16)
            v.var_desc=k.ljust(16)
            return v
        elif k in ['IONE', 'ITWO', 'KCELL','FLOW','PLMHT']:
            data_type={'IONE':'i', 'ITWO':'i', 'KCELL':'i','FLOW':'f','PLMHT':'f'}[k]
            v=self.createVariable(k,data_type,('TSTEP','NSTK'))
            v.units={'IONE':'#', 'ITWO':'#', 'KCELL':'#', 'FLOW':'m**3/hr', 'PLMHT':'m'}[k]
            v.long_name=k.ljust(16)
            v.var_desc=k.ljust(16)            
            vals = self.__hourly_stk_props[:,:,['IONE','ITWO','KCELL','FLOW','PLMHT'].index(k)]
            v[:] = vals.view('>' + data_type)
            return v
        elif k in self.__spc_names:
            v=PseudoNetCDFVariable(self,k,'f',('TSTEP','NSTK'),values=self.__emiss_data[:,self.__getspcidx(k),:])
            v.units='mole/hr'.ljust(16)
            v.long_name=k.ljust(16)
            v.var_desc=k.ljust(16)
            return v
        else:
            raise KeyError, "Unknown key %s" % k

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testPT(self):
        import PseudoNetCDF.testcase
        emissfile=point_source(PseudoNetCDF.testcase.camxfiles_paths['point_source'])
        v = emissfile.variables['NO2']
        self.assert_((v[:] == np.array([  0.00000000e+00, 3.12931000e+02, 1.23599997e+01, 0.00000000e+00, 5.27999992e+01, 0.00000000e+00, 3.12931000e+02, 1.23599997e+01, 0.00000000e+00, 5.27999992e+01], dtype = 'f').reshape(2,5)).all())

if __name__ == '__main__':
    unittest.main()
