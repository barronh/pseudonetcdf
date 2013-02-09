__all__=['lateral_boundary']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- lateral_boundary Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              lateral_boundary files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan,testing

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoIOAPIVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime
from PseudoNetCDF.camxfiles.units import get_uamiv_units

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]
class lateral_boundary(PseudoNetCDFFile):
    """
    lateral_boundary provides a PseudoNetCDF interface for CAMx
    lateral_boundary files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> lateral_boundary_path = 'camx_lateral_boundary.bin'
        >>> rows,cols = 65,83
        >>> lateral_boundaryfile = lateral_boundary(lateral_boundary_path,rows,cols)
        >>> lateral_boundaryfile.variables.keys()
        ['TFLAG', 'O3', 'NO', 'NO2', ...]
        >>> tflag = lateral_boundaryfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = lateral_boundaryfile.variables['O3']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> lateral_boundaryfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
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
        
        self.createDimension('DATE-TIME', 2)

        self.__readheader()
        
        # Add IOAPI metavariables
        nlays=self.NLAYS=len(self.dimensions['LAY'])
        nrows=self.NROWS=len(self.dimensions['ROW'])
        ncols=self.NCOLS=len(self.dimensions['COL'])
        nvars=self.NVARS=len(self.dimensions['VAR'])
        nsteps=self.NSTEPS=len(self.dimensions['TSTEP'])
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

        nx=self.__grid_hdr['nx'][0]
        ny=self.__grid_hdr['ny'][0]
        nz=max(self.__grid_hdr['nz'],array([1]))[0]

        offset+=self.__grid_hdr.dtype.itemsize*self.__grid_hdr.size
        self.__cell_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__cell_hdr_fmt,shape=1,offset=offset)
        
        
        offset+=self.__cell_hdr.dtype.itemsize*self.__cell_hdr.size+4
        self.__spc_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__spc_fmt,shape=nspec,offset=offset)
        
        offset+=self.__spc_hdr.dtype.itemsize*self.__spc_hdr.size+4
        self.__boundary_def = {}
        for bkey, bdim in [('west', ny), ('east', ny), ('south', nx), ('north', nx)]:
            __bound_fmt = dtype(dict(names=['SPAD','ione','iedge','ncell','edgedata','EPAD'],formats=['>i', '>i', '>i', '>i','(%d,4)>f' % bdim, '>i']))
            self.__boundary_def[bkey] = memmap(self.__rffile, mode = self.__mode, dtype = __bound_fmt, shape = 1, offset = offset)
            offset+=__bound_fmt.itemsize
        
        date_time_fmt=dtype(dict(names=['SPAD','BDATE','BTIME','EDATE','ETIME','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))
        date_time_block_size=6
        spc_we_fmt=dtype(dict(names=['SPAD','IONE','SPC','IEDGE', 'DATA','EPAD'],formats=['>i','>i','(10,4)>S1','>i','(%d,%d)>f' % (ny,nz),'>i']))
        spc_sn_fmt=dtype(dict(names=['SPAD','IONE','SPC','IEDGE','DATA','EPAD'],formats=['>i','>i','(10,4)>S1','>i','(%d,%d)>f' % (nx,nz),'>i']))
        
        spc_lat_fmt = dtype(dict(names = ['WEST', 'EAST', 'SOUTH', 'NORTH'], formats = [spc_we_fmt, spc_we_fmt, spc_sn_fmt, spc_sn_fmt,])) 
        # Get species names from spc_hdr
        self.__var_names__=[]
        self.__spc_names__=[''.join(spc[:,0]).strip() for spc in self.__spc_hdr]
        for spc in self.__spc_names__:
            for bkey in ['WEST_', 'EAST_', 'SOUTH_', 'NORTH_']:
                self.__var_names__.append(bkey+spc)

        spc_lat_block_size = spc_lat_fmt.itemsize / 4
        data_block_fmt=dtype(dict(names=['DATE']+self.__spc_names__,formats=[date_time_fmt]+[spc_lat_fmt]*nspec))
        
        data_block_size=date_time_block_size+nspec*spc_lat_block_size
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

    def __variables(self,k):
        spc_index=self.__var_names__.index(k)
        edgename = k.split('_')[0]
        spcname = k[len(edgename)+1:]
        dimensions=('TSTEP','LAY','ROW','COL')
        outvals=self.__memmap__[spcname][edgename]['DATA']
        units = get_uamiv_units(self.NAME, k)
        
        return PseudoIOAPIVariable(self,k,'f',dimensions,values=outvals, units = units)

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
    def testLB(self):
        import PseudoNetCDF.testcase
        latfile=lateral_boundary(PseudoNetCDF.testcase.CAMxBoundaryConditions)
        latfile.variables['TFLAG']
        aassert = testing.assert_array_almost_equal        
        aassert(latfile.variables['WEST_O3'].mean(0).mean(0), array([ 0.04808197,  0.0492915 ,  0.05001351,  0.05058479,  0.05154518,
        0.05243689,  0.05290997,  0.05321919,  0.05370211,  0.05439519,
        0.05534975,  0.05724045,  0.06059229,  0.06795996,  0.08714935,
        0.12116305],dtype='f'))
        aassert(latfile.variables['EAST_O3'].mean(0).mean(0), array([ 0.05083468,  0.05205313,  0.05279936,  0.05333328,  0.05405153,
        0.05494502,  0.05536552,  0.05566284,  0.05602615,  0.05638235,
        0.05675893,  0.05746903,  0.05893617,  0.06294519,  0.07016463,
        0.07722741],dtype='f'))
        aassert(latfile.variables['SOUTH_O3'].mean(0).mean(0), array([ 0.03825374,  0.03862067,  0.03887746,  0.03906765,  0.03936813,
        0.03965909,  0.03992134,  0.04021677,  0.04133762,  0.04258892,
        0.04378233,  0.04545989,  0.04844284,  0.05464299,  0.06876393,
        0.08342279],dtype='f'))
        aassert(latfile.variables['NORTH_O3'].mean(0).mean(0), array([ 0.03276304,  0.03307483,  0.0332762 ,  0.03341749,  0.03368971,
        0.03414582,  0.03467977,  0.03528095,  0.03693161,  0.03923673,
        0.04262673,  0.048292  ,  0.05590018,  0.0680645 ,  0.08702289,
        0.11128122],dtype='f'))
       
if __name__ == '__main__':
    lb = lateral_boundary('../../testcase/utils/CAMx/lateral_boundary/BC.vistas_2002gt2a_STL_36_68X68_16L.2002154')
    #unittest.main()