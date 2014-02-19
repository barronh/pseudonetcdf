__all__=['one3d']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- one3d Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx generic
              one 3d variable files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class one3d(PseudoNetCDFFile):
    """
    one3d provides a PseudoNetCDF interface for CAMx
    one3d files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> one3d_path = 'camx_one3d.bin'
        >>> rows,cols = 65,83
        >>> one3dfile = one3d(one3d_path,rows,cols)
        >>> one3dfile.variables.keys()
        ['TFLAG', 'UNKNOWN']
        >>> tflag = one3dfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = one3dfile.variables['UNKNOWN']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> one3dfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    

    
    id_fmt="fi"
    data_fmt="f"
    var_name="UNKNOWN"
    units="UNKNOWN"
    def __init__(self,rf,rows,cols):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        
        self.rffile=rf

        self.__memmap=memmap(self.rffile,'>f','r',offset=0)

        self.__record_items=rows*cols+4

        self.__records=self.__memmap.shape[0]/self.__record_items
        time_date=array(self.__memmap.reshape(self.__records,self.__record_items)[:,1:3])

        lays=where(time_date!=time_date[newaxis,0])[0][0]

        new_hour=slice(0,None,lays)

        dates=time_date[:,1].view('>i')
        times=time_date[:,0]

        self.__tflag=array([dates[new_hour],times[new_hour]],dtype='>f').swapaxes(0,1)
        time_steps=self.__records/lays

        self.createDimension('VAR', 1)
        self.createDimension('TSTEP', time_steps)
        self.createDimension('COL', cols)
        self.createDimension('ROW', rows)
        self.createDimension('LAY', lays)
        self.createDimension('DATE-TIME', 2)
        
        self.FTYPE=1
        self.NVARS=1
        self.NCOLS=cols
        self.NROWS=rows
        self.NLAYS=lays
        self.NTHIK=1
        
        self.variables=PseudoNetCDFVariables(self.__variables,[self.var_name])
        v=self.variables['TFLAG']=ConvertCAMxTime(self.__tflag[:,0],self.__tflag[:,1],1)
        self.SDATE,self.STIME=v[0,0,:]

    def __decorator(self,name,pncfv):
        decor=lambda *args: dict(units=self.units, var_desc=self.var_name.ljust(16), long_name=self.var_name.ljust(16))
        for k,v in decor(name).iteritems():
            setattr(pncfv,k,v)
        return pncfv
        
    def __variables(self,k):
        tsteps=len(self.dimensions['TSTEP'])
        lays=len(self.dimensions['LAY'])
        rows=len(self.dimensions['ROW'])
        cols=len(self.dimensions['COL'])
        return self.__decorator(k,PseudoNetCDFVariable(self,k,'f',('TSTEP','LAY','ROW','COL'),values=self.__memmap.reshape(self.__records,self.__record_items)[:,3:-1].reshape(tsteps,lays,rows,cols)))

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testKV(self):
        import PseudoNetCDF.testcase
        vdfile=one3d(PseudoNetCDF.testcase.camxfiles_paths['vertical_diffusivity'],4,5)
        vdfile.variables['TFLAG']
        self.assert_((vdfile.variables['UNKNOWN']==array([1.00000000e+00, 4.76359320e+00, 1.92715893e+01, 1.52158489e+01, 7.20601225e+00, 1.84097159e+00, 2.63084507e+01, 1.27621298e+01, 1.06348248e+01, 2.22587357e+01, 1.00000000e+00, 1.69009724e+01, 1.00000000e+00, 2.23075104e+01, 1.27485418e+01, 1.45508013e+01, 1.45637455e+01, 2.95294094e+01, 3.02676849e+01, 2.84974957e+01, 1.00000000e+00, 6.93706131e+00, 3.07418957e+01, 3.41300621e+01, 1.66266994e+01, 1.00000000e+00, 2.53920174e+01, 1.92539787e+01, 2.32906532e+01, 5.96702042e+01, 1.00000000e+00, 2.24458847e+01, 1.00000000e+00, 5.45038452e+01, 3.45825729e+01, 3.43578224e+01, 3.76071548e+01, 6.76799850e+01, 7.33648529e+01, 7.35801239e+01, 1.00000000e+00, 6.38178444e+00, 4.39327278e+01, 5.00754166e+01, 2.44474106e+01, 1.00000000e+00, 3.51700935e+01, 2.71137428e+01, 3.40312347e+01, 8.85775909e+01, 1.00000000e+00, 3.13994522e+01, 1.00000000e+00, 8.07169266e+01, 5.12892876e+01, 5.05734329e+01, 5.56603966e+01, 1.00394188e+02, 1.08980370e+02, 1.09251083e+02, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.64916098e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.22205174e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.05408611e+01, 1.26687222e+01, 1.21386652e+01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.93040049e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.39350021e+00, 1.02697349e+00, 1.00000000e+00, 1.00000000e+00, 1.82250175e+01, 2.90407104e+01, 2.83827496e+01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 2.02609706e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.44422662e+00, 1.02998519e+00, 1.00000000e+00, 1.00000000e+00, 2.60322971e+01, 4.26534195e+01, 4.17046585e+01],dtype='f').reshape(2,3,4,5)).all())
    
       
if __name__ == '__main__':
    unittest.main()
