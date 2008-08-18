HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['point_source']
#Distribution packages
import unittest
import struct
from warnings import warn
#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd
from pyPA.utils.FortranFileUtil import OpenRecordFile,Int2Asc
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class point_source(PseudoNetCDFFile):
    """
    This class is intended to provide an interface to the
    point source emission CAMx format.
    """
    __emiss_hdr_fmt=dtype(dict(names=['SPAD','name','note','ione','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>10i','>60i','>i','>i','>i','>f','>i','>f','>i']))
    __grid_hdr_fmt=dtype(dict(names=['SPAD','rdum1','rdum2','iutm','xorg','yorg','delx','dely','nx','ny','nz', 'idum1','idum2','rdum3','rdum4','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))
    __cell_hdr_fmt=dtype(dict(names=['SPAD','ione1','ione2','nx','ny','EPAD'],formats=['>i','>i','>i','>i','>i','>i']))
    
    __spc_hdr_fmt=dtype('>10i')
    __time_hdr_fmt=dtype(dict(names=['SPAD','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))
    __nstk_hdr_fmt=dtype(dict(names=['SPAD','ione','nstk','EPAD'],formats=['>i','>i','>i','>i']))
    
    
    __stk_prop_fmt=dtype(dict(names=['XSTK','YSTK','HSTK','DSTK','TSTK','VSTK'],formats=['>f','>f','>f','>f','>f','>f']))
    
    __stk_time_prop_fmt=dtype(dict(names=['idum1','idum2','KCELL','FLOW','PLMHT'],formats=['>i','>i','>f','>f','>f']))

    def __init__(self,rf):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.__rffile=rf
        
        
        self.dimensions={}
        self.variables={}
        self.__memmap=memmap(self.__rffile,'>f','r')
        self.__globalheader()
        varkeys=('ETFLAG','TFLAG','XSTK','YSTK','HSTK','DSTK','TSTK','VSTK','KCELL','FLOW','PLMHT','NSTKS')+tuple([i.strip() for i in self.__spc_names])
        self.variables=PseudoNetCDFVariables(self.__variables,varkeys)
        self.__time_stks()
        self.createDimension('STK',self.__nstk_hdr['nstk'])
        
        
    
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

        self.__cell_hdr=self.__memmap[offset:offset+self.__cell_hdr_fmt.itemsize/4].view(self.__cell_hdr_fmt)
        offset+=self.__cell_hdr.nbytes/4+1

        nspec=self.__emiss_hdr['nspec'][0]


        self.__spc_hdr=self.__memmap[offset:offset+self.__spc_hdr_fmt.itemsize/4*nspec].reshape(nspec,self.__spc_hdr_fmt.itemsize/4).view('>i')
        offset+=self.__spc_hdr.nbytes/4+1

        self.__spc_names=[Int2Asc(spc).strip() for spc in self.__spc_hdr]
        self.__nstk_hdr=self.__memmap[offset:offset+self.__nstk_hdr_fmt.itemsize/4].view(self.__nstk_hdr_fmt)
        offset+=self.__nstk_hdr.nbytes/4+1

        self.dimensions['NSTK']=nstk=self.__nstk_hdr['nstk']
        

        self.__nstk_hdr['nstk']
        self.__stk_props=self.__memmap[offset:offset+self.__stk_prop_fmt.itemsize/4*nstk].reshape(nstk,self.__stk_prop_fmt.itemsize/4).view(self.__stk_prop_fmt)
        offset+=self.__stk_props.nbytes/4+1

        self.__data_start=offset
        
        self.SDATE=self.__emiss_hdr['ibdate']
        self.STIME=self.__emiss_hdr['btime']
        self.EDATE=self.__emiss_hdr['iedate']
        self.ETIME=self.__emiss_hdr['etime']
        
        self.TSTEP=timediff((self.SDATE,self.STIME),(self.EDATE,self.ETIME))

    def __getspcidx(self,spc):
        return self.__spc_names.index(spc)

    def __time_stks(self):
        i=offset=0
        nspcs=len(self.__spc_names)
        nstks=self.dimensions['NSTK']
        date_block_size=6
        stk_block_size=4
        stk_props_size=2+nstks*5
        emiss_block_size=nspcs*(nstks+13)
        hour_block_size=date_block_size+stk_block_size+stk_props_size+emiss_block_size
        data=self.__memmap[self.__data_start:]
        data=data.reshape(data.size/hour_block_size,hour_block_size)
        ntimes=data.shape[0]
        self.createDimension('TSTEP',ntimes)

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
        
        self.variables['TFLAG']=ConvertCAMxTime(bdates,btimes,self.dimensions['VAR'])
        self.variables['ETFLAG']=ConvertCAMxTime(edates,etimes,self.dimensions['VAR'])
        v=self.variables['NSTKS']=PseudoNetCDFVariable(self,'NSTKS','i',('TSTEP',),array(nstk_hdr))
        v.units='#'.ljust(16)
        v.long_name='NSTKS'.ljust(16)
        v.var_desc=v.long_name
    
    def __variables(self,k):
        if k in ['TFLAG','ETFLAG','NSTKS']:
            return self.variables[k]
        elif k in ['XSTK','YSTK','HSTK','DSTK','TSTK','VSTK']:
            v=PseudoNetCDFVariable(self,k,'f',('NSTK',),self.__stk_props[k])
            v.units={'XSTK':'m','YSTK':'m','HSTK':'m','DSTK':'m','TSTK':'K','VSTK':'m/h'}[k]
            v.long_name=k.ljust(16)
            v.var_desc=k.ljust(16)
            return v
        elif k in ['KCELL','FLOW','PLMHT']:
            data_type={'KCELL':'i','FLOW':'f','PLMHT':'f'}[k]
            v=self.createVariable(k,data_type,('TSTEP','NSTK'))
            v.units={'KCELL':'#','FLOW':'m**3/hr','PLMHT':'m'}[k]
            v.long_name=k.ljust(16)
            v.var_desc=k.ljust(16)            
            v.assignValue(self.__hourly_stk_props[:,:,['',' ','KCELL','FLOW','PLMHT'].index(k)])
            return v
        elif k in self.__spc_names:
            v=PseudoNetCDFVariable(self,k,'f',('TSTEP','NSTK'),self.__emiss_data[:,self.__getspcidx(k),:])
            v.units='mole/hr'.ljust(16)
            v.long_name=k.ljust(16)
            v.var_desc=k.ljust(16)
            return v
        else:
            raise KeyError, "Unknown key %s" % k

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testPT(self):
        emissfile=point_source('../../../../testdata/ei/camx_cb4_ei_el.20000825.hgb8h.base1b.psito2n2')
        for k in emissfile.variables.keys():
            v=emissfile.variables[k]
            if k in ['TFLAG','ETFLAG']:
                print k,v[[0,-1],0,:]
            else:
                print k,v.min(),v.mean(),v.max()
        import pdb; pdb.set_trace()
        self.assert_(1==2)

if __name__ == '__main__':
    unittest.main()
