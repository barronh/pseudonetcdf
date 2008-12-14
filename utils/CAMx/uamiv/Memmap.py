HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['uamiv']
#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd
from pyPA.utils.FortranFileUtil import OpenRecordFile
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class uamiv(PseudoNetCDFFile):
    """
    This class is intended to provide an interface to the
    low level (gridded emission) or initial conditions 
    inputs and instantaneous or average outputs of the 
    CAMx model
    """
    
    __emiss_hdr_fmt=dtype(dict(names=['SPAD','name','note','ione','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','(10,4)>S1','(60,4)>S1','>i','>i','>i','>f','>i','>f','>i']))
    __grid_hdr_fmt=dtype(dict(names=['SPAD','rdum1','rdum2','iutm','xorg','yorg','delx','dely','nx','ny','nz','idum1','idum2','rdum3','rdum4','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))
    __cell_hdr_fmt=dtype(dict(names=['SPAD','ione1','ione2','nx','ny','EPAD'],formats=['>i','>i','>i','>i','>i','>i']))
    __time_hdr_fmt=dtype(dict(names=['SPAD','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))
    __spc_fmt=dtype("(10,4)>S1")
    
    __ione=1
    __idum=0
    __rdum=0.
    def __init__(self,rf,mode='r'):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.__rffile=rf
        self.__mode=mode
        # Establish dimensions
        self.dimensions={'DATE-TIME': 2 }

        self.__readheader()
        
        # Add IOAPI metavariables
        nlays=self.NLAYS=self.dimensions['LAY']
        nrows=self.NROWS=self.dimensions['ROW']
        ncols=self.NCOLS=self.dimensions['COL']
        nvars=self.NVARS=self.dimensions['VAR']
        nsteps=self.NSTEPS=self.dimensions['TSTEP']
        setattr(self,'VAR-LIST',"".join([i.ljust(16) for i in self.__var_names__]+['TFLAG'.ljust(16)]))
        self.GDTYP=2
        self.NAME="".join(self.__emiss_hdr['name'][0,:,0])

        # Create variables
        self.variables=PseudoNetCDFVariables(self.__variables,self.__var_names__+['TFLAG','ETFLAG'])
        self.variables['TFLAG']=ConvertCAMxTime(self.__memmap__['DATE']['BDATE'],self.__memmap__['DATE']['BTIME'],self.NVARS)
        self.variables['ETFLAG']=ConvertCAMxTime(self.__memmap__['DATE']['EDATE'],self.__memmap__['DATE']['BTIME'],self.NVARS)
        
        self.SDATE,self.STIME=self.variables['TFLAG'][0,0,:]
        
    def __checkfilelen(self):
        f=file(self.__rffile,'rb')
        f.seek(0,2)
        flen=f.tell()
        f.close()
        return flen

    def __readheader(self):
        offset=0
        self.__emiss_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__emiss_hdr_fmt,shape=1,offset=offset)
        nspec=self.__emiss_hdr['nspec'][0]
        offset+=self.__emiss_hdr.dtype.itemsize*self.__emiss_hdr.size

        self.__grid_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__grid_hdr_fmt,shape=1,offset=offset)

        self.XORIG=self.__grid_hdr['xorg'][0]
        self.YORIG=self.__grid_hdr['yorg'][0]
        self.XCELL=self.__grid_hdr['delx'][0]
        self.YCELL=self.__grid_hdr['dely'][0]

        nx=self.__grid_hdr['nx']
        ny=self.__grid_hdr['ny']
        nz=max(self.__grid_hdr['nz'],1)

        offset+=self.__grid_hdr.dtype.itemsize*self.__grid_hdr.size
        self.__cell_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__cell_hdr_fmt,shape=1,offset=offset)
        
        
        offset+=self.__cell_hdr.dtype.itemsize*self.__cell_hdr.size+4
        self.__spc_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__spc_fmt,shape=nspec,offset=offset)
        
        offset+=self.__spc_hdr.dtype.itemsize*self.__spc_hdr.size+4
        
        
        date_time_fmt=dtype(dict(names=['SPAD','BDATE','BTIME','EDATE','ETIME','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))
        date_time_block_size=6
        spc_1_lay_fmt=dtype(dict(names=['SPAD','IONE','SPC','DATA','EPAD'],formats=['>i','>i','(10,4)>S1','(%d,%d)>f' % (ny,nx),'>i']))
        spc_1_lay_block_size=13+nx*ny
        spc_3d_fmt=dtype((spc_1_lay_fmt,(nz,)))
        spc_3d_block_size=spc_1_lay_block_size*nz
        
        # Get species names from spc_hdr
        self.__var_names__=[''.join(spc[:,0]).strip() for spc in self.__spc_hdr]

        data_block_fmt=dtype(dict(names=['DATE']+self.__var_names__,formats=[date_time_fmt]+[spc_3d_fmt]*nspec))
        
        data_block_size=date_time_block_size+nspec*nz*spc_1_lay_block_size
        f=file(self.__rffile)
        f.seek(0,2)
        size=f.tell()
        f.close()
        del f
        ntimes=float(size-offset)/4./data_block_size
        if int(ntimes)!=ntimes:
            raise ValueError, "Not an even number of times %f" % ntimes
        ntimes=int(ntimes)
        
        self.createDimension('LAY',nz)
        self.createDimension('COL',nx)
        self.createDimension('ROW',ny)
        self.createDimension('TSTEP',ntimes)
        self.createDimension('VAR',nspec)
        
        self.__memmap__=memmap(self.__rffile,mode=self.__mode,dtype=data_block_fmt,offset=offset)

    def __decorator(self,name,pncfv):
        if name=='EMISSIONS ':
            decor=lambda spc: dict(units='umol/hr', var_desc=spc.ljust(16), long_name=spc.ljust(16))
        else:
            decor=lambda spc: dict(units='ppm', var_desc=spc.ljust(16), long_name=spc.ljust(16))

        for k,v in decor(name).iteritems():
            setattr(pncfv,k,v)
        return pncfv
        
    def __variables(self,k):
        spc_index=self.__var_names__.index(k)
        dimensions=('TSTEP','LAY','ROW','COL')
        outvals=self.__memmap__[k]['DATA']
        return self.__decorator(k,PseudoNetCDFVariable(self,k,'f',dimensions,values=outvals))

    def sync(self):
        pass
    
    def close(self):
        self.sync()
        self.__memmap__.close()
        

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
    def testGE(self):
        import pyPA.testcase
        emissfile=uamiv(pyPA.testcase.CAMxAreaEmissions)
        emissfile.variables['TFLAG']
        v=emissfile.variables['NO']
        self.assert_((emissfile.variables['NO'].mean(1).mean(1).mean(1)==array([  52.05988312,   51.58646774,   51.28796387,   55.63090134,
         63.95315933,  105.3456192 ,  158.26776123,  152.04057312,
         147.32403564,  154.80661011,  164.03274536,  171.88658142,
         174.36567688,  180.03359985,  173.81938171,  180.50257874,
         178.56637573,  161.35736084,  110.38669586,   97.90225983,
         89.08138275,   81.10474396,   73.36611938,   58.82622528],dtype='f')).all())

    def testAvg(self):
        import pyPA.testcase
        emissfile=uamiv(pyPA.testcase.CAMxAverage)
        emissfile.variables['TFLAG']
        v=emissfile.variables['NO']
        self.assert_((v.mean(1).mean(1).mean(1)==array([  9.44490694e-06,   2.17493564e-07,   6.08432686e-07,   9.48155161e-07,
         1.15099192e-05,   1.02132122e-04,   2.57815613e-04,   3.35910037e-04,
         3.17813188e-04,   2.51695659e-04,   1.85225872e-04,   1.40698961e-04,
         1.16110547e-04,   1.04519037e-04,   1.00367179e-04,   9.81789271e-05,
         8.98482831e-05,   6.31201983e-05,   2.18762198e-05,   1.78832056e-06,
         1.20556749e-07,   1.57714638e-07,   1.82648236e-07,   2.02759026e-07],dtype='f')).all())

    def testInst(self):
        from warnings import warn
        warn("Instantaneous file test not implemented")
       
if __name__ == '__main__':
    unittest.main()
