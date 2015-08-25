__all__=['finst']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- finst Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx fine
              instantaneus/avg files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

from warnings import warn
#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan,fromfile

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd
from PseudoNetCDF.camxfiles.units import get_uamiv_units
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class finst(PseudoNetCDFFile):
    """
    finst provides a PseudoNetCDF interface for CAMx
    finst files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> finst_path = 'camx_finst.bin'
        >>> finstfile = finst(finst_path)
        >>> finstfile.variables.keys()
        ['TFLAG', 'O3', 'NO', 'NO2', ...]
        >>> v = finstfile.variables['V']
        >>> tflag = finstfile.variables['TFLAG']
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
        >>> finstfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    
    __messagefmt='>240S'
    __nestspcfmt=dtype(dict(names=['nnest','nspec'],formats=['>i','>i']))
    __nestparmsfmt=dtype(dict(names=['SPAD','ibeg','jbeg','iend','jend','mesh','ione','nx','ny','nz','iparent','ilevel','EPAD'],formats=['>i']*13))
    __time_hdr_fmt=dtype(dict(names=['time','date'],formats=['>f','>i']))
    __spc_fmt=dtype(">10S")
    __ione=1
    __idum=0
    __rdum=0.
    __buffersize=4
    def __init__(self,rf,mode='r'):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.__rffile=rf

        # Establish dimensions
        self.createDimension('DATE-TIME', 2)

        self.__readheader(mode)
        
        # Add IOAPI metavariables
        nlays=self.NLAYS=len(self.dimensions['LAY'])
        nrows=self.NROWS=len(self.dimensions['ROW'])
        ncols=self.NCOLS=len(self.dimensions['COL'])
        nvars=self.NVARS=len(self.dimensions['VAR'])
        nsteps=self.NSTEPS=len(self.dimensions['TSTEP'])
        setattr(self,'VAR-LIST',"".join([i.ljust(16) for i in self.__var_names__] + ['TFLAG'.ljust(16)]))
        self.GDTYP=2

        # Create variables
        self.variables=PseudoNetCDFVariables(self.__variables,self.__var_names__+['TFLAG'])

        # Initialize time maps
        date=self.__memmap__['DATE']['DATE']
        time=self.__memmap__['DATE']['TIME']
        self.variables['TFLAG']=ConvertCAMxTime(date,time,self.NVARS)
        self.SDATE,self.STIME=self.variables['TFLAG'][0,0,:]
        

    def __checkfilelen(self):
        f=open(self.__rffile,'rb')
        f.seek(0,2)
        flen=f.tell()
        f.close()
        return flen

    def __readheader(self,mode):
        start=0
        end=0
        #start+=self.__buffersize
        #end=self.__messagefmt.itemsize/4
        #self.MESSAGE=self.__memmap__[start:end].view(self.__messagefmt)
        f=open(self.__rffile)
        f.seek(92)
        self.NNEST,self.NSPEC=fromfile(f,'>i',2)
        f.seek(8,1)

        self.SPECIES=fromfile(f,self.__spc_fmt,self.NSPEC)
        
        offset=f.tell()+4
        f.close()
        del f
        self.__var_names__=[i.strip() for i in self.SPECIES.tolist()]

        self.NEST_PARMS=memmap(self.__rffile,dtype=self.__nestparmsfmt,mode=mode,offset=offset,shape=self.NNEST)
        
        offset+=self.__nestparmsfmt.itemsize*self.NNEST
        
        date_fmt=dtype(dict(names=['SPAD','TIME','DATE','EPAD'],formats=['>i','>f','>i','>i']))
        spc_1_lay_fmt=[]
        for i in range(self.NNEST):
            spc_1_lay_fmt.append( dtype( dict( names = ['SPAD', 'DATA', 'EPAD'], formats = ['>i', '(%d,%d)>f' % (self.NEST_PARMS[i]['ny'], self.NEST_PARMS[i]['nx']), '>i'] ) ) )
        
        grid_fmt=[]
        for i in range(self.NNEST):
            grid_fmt.append( dtype( dict( names=self.__var_names__, formats=[dtype( ( spc_1_lay_fmt[i], ( int(self.NEST_PARMS[i]['nz']), ) ) )]*self.NSPEC ) ) )
        
        
        self.__memmap__=memmap(self.__rffile, mode=mode, offset=offset, dtype=dtype(dict(names=['DATE']+['grid%d' % i for i in range(self.NNEST)],formats=[date_fmt]+grid_fmt)))
        
        ntimes=self.__memmap__.shape[0]
        if int(ntimes)!=ntimes:
            raise ValueError, "Not an even number of times"

        self.createDimension('TSTEP',ntimes)
        self.createDimension('VAR',self.NSPEC)
        self.chooseGrid(0)
        
        
        
    def chooseGrid(self,ngrid):
        self.createDimension('LAY',self.NEST_PARMS[ngrid]['nz'])
        self.createDimension('COL',self.NEST_PARMS[ngrid]['nx'])
        self.createDimension('ROW',self.NEST_PARMS[ngrid]['ny'])
        self.CURRENT_GRID='grid%d' % ngrid

    def sync(self):
        pass
    
    def close(self):
        self.sync()
        self.__memmap__.close()
        
    def __variables(self,k):
        dimensions=('TSTEP','LAY','ROW','COL')
        ntimes=len(self.dimensions['TSTEP'])
        nx=len(self.dimensions['COL'])
        ny=len(self.dimensions['ROW'])
        nz=len(self.dimensions['LAY'])
        
        outvals=self.__memmap__[self.CURRENT_GRID][k]['DATA'][:,:,:,:]
        unit = get_uamiv_units('INSTANT   ', k)
        return PseudoNetCDFVariable(self,k,'f',dimensions,values=outvals,units=unit)

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testFinst(self):
        import PseudoNetCDF.testcase
        finstfile=finst(PseudoNetCDF.testcase.camxfiles_paths['finst'])
        vars=[i for i in finstfile.variables.keys() if i!='TFLAG']
        for var in vars:
            v=finstfile.variables[var]
            warn("Test case is not fully implemented.  Review values for 'reasonability.'")
            datetime=finstfile.variables['TFLAG']
            date=datetime[:,0,0]
            time=datetime[:,0,1]
            minv=v.min(1).min(1).min(1)
            meanv=v.mean(1).mean(1).mean(1)
            maxv=v.max(1).max(1).max(1)
            for i,(d,t) in enumerate(datetime[:,0,:]):
                print var,d,t,minv[i],meanv[i],maxv[i]
if __name__ == '__main__':
    unittest.main()
