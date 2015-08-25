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
    

    __emiss_hdr_fmt=dtype(dict(names=['SPAD','name','note','itzon','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','(10,4)>S1','(60,4)>S1','>i','>i','>i','>f','>i','>f','>i']))
    __grid_hdr_fmt=dtype(dict(names=['SPAD','plon','plat','iutm','xorg','yorg','delx','dely','nx','ny','nz','iproj','istag','tlat1','tlat2','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))
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
        self.NOTE="".join(self.__emiss_hdr['note'][0,:,0])
        self.ITZON=self.__emiss_hdr['itzon'][0]

        # Create variables
        self.variables=PseudoNetCDFVariables(self.__variables,self.__var_names__+['TFLAG','ETFLAG'])
        self.variables['TFLAG']=ConvertCAMxTime(self.__memmap__['DATE']['BDATE'],self.__memmap__['DATE']['BTIME'],self.NVARS)
        self.variables['ETFLAG']=ConvertCAMxTime(self.__memmap__['DATE']['EDATE'],self.__memmap__['DATE']['BTIME'],self.NVARS)
        
        self.SDATE,self.STIME=self.variables['TFLAG'][0,0,:]
        
    def __checkfilelen(self):
        f=open(self.__rffile,'rb')
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
        self.PLON = plon = self.__grid_hdr['plon'][0]
        self.PLAT = plat = self.__grid_hdr['plat'][0]
        self.TLAT1 = tlat1 = self.__grid_hdr['tlat1'][0]
        self.TLAT2 = tlat2 = self.__grid_hdr['tlat2'][0]
        self.IUTM = iutm = self.__grid_hdr['iutm'][0]
        self.ISTAG = istag = self.__grid_hdr['istag'][0]
        self.CPROJ = cproj = self.__grid_hdr['iproj'][0]

        # Map CAMx projection constants to IOAPI
        GDTYPE = self.GDTYP={0: 1, 1: 5, 2: 2, 3: 6}[self.__grid_hdr['iproj'][0]]
        self.XCENT = self.__grid_hdr['plon'][0]
        self.YCENT = self.__grid_hdr['plat'][0]
        if GDTYPE in (1, 2):
            self.P_ALP = self.__grid_hdr['tlat1'][0]
            self.P_BET = self.__grid_hdr['tlat2'][0]
            self.P_GAM = self.__grid_hdr['plon'][0]
        elif GDTYPE == 5:
            self.P_ALP = self.__grid_hdr['iutm'][0]
            self.P_BET = 0.
            self.P_GAM = 0.
        elif GDTYPE == 6:
            self.P_ALP = {90: 1, -90: -1}[self.__grid_hdr['plat'][0]]
            self.P_BET = self.__grid_hdr['tlat1'][0]
            self.P_GAM = self.__grid_hdr['plon'][0]
        else:
            raise ValueError('Unknown projection')            

        nx=self.__grid_hdr['nx'][0]
        ny=self.__grid_hdr['ny'][0]
        nz=max(self.__grid_hdr['nz'],array([1]))[0]

        offset+=self.__grid_hdr.dtype.itemsize*self.__grid_hdr.size
        self.__cell_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__cell_hdr_fmt,shape=1,offset=offset)
        
        
        offset+=self.__cell_hdr.dtype.itemsize*self.__cell_hdr.size+4
        self.__spc_hdr=memmap(self.__rffile,mode=self.__mode,dtype=self.__spc_fmt,shape=nspec,offset=offset)
        
        offset+=self.__spc_hdr.dtype.itemsize*self.__spc_hdr.size+4
        self._boundary_def = {}
        for bkey, bdim in [('WEST', ny), ('EAST', ny), ('SOUTH', nx), ('NORTH', nx)]:
            __bound_fmt = dtype(dict(names=['SPAD','ione','iedge','ncell','edgedata','EPAD'],formats=['>i', '>i', '>i', '>i','(%d,%d)>i' % (bdim, 4), '>i']))
            self._boundary_def[bkey] = memmap(self.__rffile, mode = self.__mode, dtype = __bound_fmt, shape = 1, offset = offset)
            assert(self._boundary_def[bkey]['SPAD'] == (__bound_fmt.itemsize - 8))
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
        f=open(self.__rffile)
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
        self.createDimension('VAR',nspec*4)
        
        self.__memmap__=memmap(self.__rffile,mode=self.__mode,dtype=data_block_fmt,offset=offset)

    def __variables(self,k):
        spc_index=self.__var_names__.index(k)
        edgename = k.split('_')[0]
        spcname = k[len(edgename)+1:]
        if edgename in ('WEST', 'EAST'):
            dimensions=('TSTEP', 'ROW', 'LAY')
        else:
            dimensions=('TSTEP', 'COL', 'LAY')
        outvals=self.__memmap__[spcname][edgename]['DATA']
        units = get_uamiv_units(self.NAME, spcname)
        
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
        latfile=lateral_boundary(PseudoNetCDF.testcase.camxfiles_paths['lateral_boundary'])
        latfile.variables['TFLAG']
        aassert = testing.assert_array_almost_equal        
        aassert(latfile.variables['WEST_O3'], array([5.09086549e-02, 5.15854955e-02, 5.20327762e-02, 4.51020785e-02, 4.81765866e-02, 5.01902141e-02, 5.11084236e-02, 5.30455858e-02, 5.39216697e-02, 5.49603552e-02, 5.52441478e-02, 5.53695373e-02, 4.70436066e-02, 4.88652885e-02, 5.00564687e-02, 4.30268869e-02, 4.65402082e-02, 4.93281633e-02, 4.84376177e-02, 5.15173525e-02, 5.32316864e-02, 5.27723655e-02, 5.43252379e-02, 5.51613718e-02],dtype='f').reshape(2,4,3))
        aassert(latfile.variables['EAST_O3'], array([4.32289392e-02, 4.32576202e-02, 4.32824045e-02, 4.28533703e-02, 4.29217815e-02, 4.29747701e-02, 4.20787595e-02, 4.21368703e-02, 4.21847440e-02, 4.20016758e-02, 4.20646742e-02, 4.21276242e-02, 4.45125513e-02, 4.45423163e-02, 4.45689932e-02, 4.34477255e-02, 4.35189568e-02, 4.35763896e-02, 4.23655175e-02, 4.24236283e-02, 4.24756072e-02, 4.21473905e-02, 4.22183946e-02, 4.22968678e-02],dtype='f').reshape(2,4,3))
        aassert(latfile.variables['SOUTH_O3'], array([5.00756986e-02, 5.07836118e-02, 5.12398519e-02, 4.53975834e-02, 4.85784635e-02, 5.06106876e-02, 4.88935970e-02, 5.18746264e-02, 5.34388497e-02, 5.52802160e-02, 5.66427484e-02, 5.73641062e-02, 5.76447174e-02, 5.91668785e-02, 5.98937273e-02, 4.64089997e-02, 4.85925563e-02, 4.98120859e-02, 4.31438088e-02, 4.85876575e-02, 5.05158082e-02, 4.63983566e-02, 5.09891734e-02, 5.34970984e-02, 5.16217351e-02, 5.49310222e-02, 5.70215993e-02, 5.37858233e-02, 5.71687669e-02, 5.90757653e-02],dtype='f').reshape(2, 5, 3))
        aassert(latfile.variables['NORTH_O3'], array([4.67120931e-02, 4.72769439e-02, 4.75757420e-02, 4.63561714e-02, 4.69456315e-02, 4.72794324e-02, 4.45141681e-02, 4.53266129e-02, 4.58004810e-02, 4.20911871e-02, 4.29873653e-02, 4.35210876e-02, 3.29810269e-02, 3.38455141e-02, 3.43711227e-02, 4.65101935e-02, 4.71580848e-02, 4.75305654e-02, 4.55248095e-02, 4.63739596e-02, 4.68796641e-02, 4.44590077e-02, 4.58386838e-02, 4.66508269e-02, 4.33906466e-02, 4.47074845e-02, 4.51488644e-02, 3.23279053e-02, 3.31198536e-02, 3.36514898e-02],dtype='f').reshape(2, 5, 3))
       
if __name__ == '__main__':
    lb = lateral_boundary('../../testcase/utils/CAMx/lateral_boundary/BC.vistas_2002gt2a_STL_36_68X68_16L.2002154')
    #unittest.main()
