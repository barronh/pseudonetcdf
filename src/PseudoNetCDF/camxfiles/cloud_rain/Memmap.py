__all__ = ['cloud_rain']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- cloud_rain Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              cloud/rain files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

HeadURL = "$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
import unittest
import struct
from warnings import warn

#Site-Packages
from numpy import zeros, array, where,memmap, newaxis, dtype,nan

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff, timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan = struct.unpack('>f', '\xff\xc0\x00\x00')[0]
checkarray = zeros((1, ), 'f')
checkarray[0] = listnan
array_nan = checkarray[0]
    
class cloud_rain(PseudoNetCDFFile):
    """
    cloud_rain provides a PseudoNetCDF interface for CAMx
    cloud_rain files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> cloud_rain_path = 'cloud_rain.bin'
        >>> rows, cols = 65, 83
        >>> cloud_rainfile = cloud_rain(cloud_rain_path, rows, cols)
        >>> cloud_rainfile.variables.keys()
        ['CLOUD', 'RAIN', 'SNOW', 'GRAUPEL', 'COD', 'TFLAG']
        >>> v = cloud_rainfile.variables['CLOUD']
        >>> tflag = cloud_rainfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0, 0, :]
        array([2005185,       0])
        >>> tflag[-1, 0, :]
        array([2005185,  240000])
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> cloud_rainfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    def __init__(self, rf, rows = None,cols = None):
        f = open(rf, 'rb')
        f.seek(0, 2)
        flen = f.tell()
        offset = struct.unpack('>i', open(rf, 'r').read(4))[0] + 8
        self.__memmap = memmap(rf, '>f', 'r', offset = offset)
        ncols, nrows, nlays = struct.unpack({35:'>i15ciiii', 40:'>i20ciiii'}[offset], open(rf, 'r').read(offset))[-4:-1]
        self.createDimension('COL', ncols)
        self.createDimension('ROW', nrows)
        self.createDimension('LAY', nlays)
        header = struct.unpack({35:'>i15ciiiiifi', 40:'>i20ciiiiifi'}[offset], open(rf, 'r').read(offset + 12))
        self.FILEDESC = ''.join(header[1:1+{35: 15, 40: 20}[offset]])
        self.STIME, self.SDATE = header[-2:]
        if self.SDATE < 10000:
            self.SDATE+=2000000
        if (ncols!=cols and cols!=None) or (rows!=rows and rows!=None):
            warn('Files says cols = %d, rows = %d, and lays = %d; you said cols = %d and rows = %d' % (ncols, nrows, nlays, cols, rows))
            
        self.createDimension('DATE-TIME', 2)
        self.VERSION, varkeys = {35:('<4.3', ['CLOUD', 'PRECIP', 'COD', 'TFLAG']), 40:('4.3', ['CLOUD', 'RAIN', 'SNOW', 'GRAUPEL', 'COD', 'TFLAG'])}[offset]
        self.createDimension('TSTEP', int((flen - offset) / ((len(varkeys) - 1) * nlays * (nrows * ncols + 2) * 4 + 16)))
        self.createDimension('VAR', len(varkeys) - 1)
        
        self.NVARS = len(self.dimensions['VAR'])
        self.NLAYS = len(self.dimensions['LAY'])
        self.NROWS = len(self.dimensions['ROW'])
        self.NCOLS = len(self.dimensions['COL'])
        self.FTYPE = 1
        
        self.variables = PseudoNetCDFVariables(self.__var_get, varkeys)
        
        self.SDATE,self.STIME = self.variables['TFLAG'][0, 0, :]

    def __set_var(self, key, vals_idx):
        times = len(self.dimensions['TSTEP'])
        lays = len(self.dimensions['LAY'])
        rows = len(self.dimensions['ROW'])
        cols = len(self.dimensions['COL'])
        v = PseudoNetCDFVariable(self, key, 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), values = self.__memmap[vals_idx].reshape(times, lays, rows, cols))
        v.units = {'COD':'None'}.get(key, 'g/m**3')
        v.long_name = key
        v.var_desc = key
        self.variables[key] = v
        
    def __var_get(self, key):
        times = len(self.dimensions['TSTEP'])
        rows = len(self.dimensions['ROW'])
        cols = len(self.dimensions['COL'])
        lays = len(self.dimensions['LAY'])
        vars = len(self.variables.keys()) - 1
        hour = 1
        date = 2
        cloud = 3
        rain = 4
        snow = 5
        graupel = 6
        cod = 7
        stagger = 8
        out_idx = zeros(self.__memmap.shape,dtype = 'b')
        out_idx.reshape(times, lays * vars * (rows * cols + 2) + 4)[:,1] = hour
        out_idx.reshape(times, lays * vars * (rows * cols + 2) + 4)[:,2] = date
        
        self.variables['TFLAG'] = ConvertCAMxTime(self.__memmap[out_idx==date].view('>i'), self.__memmap[out_idx==hour], len(self.dimensions['VAR']))
        
        val_shape = out_idx.reshape(times, lays * vars * (rows * cols + 2) + 4)[:,4:].reshape(times, lays, vars, rows * cols + 2)[:, :,:,1:-1].reshape(times, lays, vars, rows, cols)
        if self.VERSION=='<4.3':
            val_shape[:, :,0, :, :] = cloud
            val_shape[:, :,1, :, :] = rain
            val_shape[:, :,2, :, :] = cod
            self.__set_var('CLOUD', out_idx==cloud)
            self.__set_var('PRECIP', out_idx==rain)
            self.__set_var('COD', out_idx==cod)
        else:
            val_shape[:, :,0, :, :] = cloud
            val_shape[:, :,1, :, :] = rain
            val_shape[:, :,2, :, :] = snow
            val_shape[:, :,3, :, :] = graupel
            val_shape[:, :,4, :, :] = cod
            self.__set_var('CLOUD', out_idx==cloud)
            self.__set_var('RAIN', out_idx==rain)
            self.__set_var('SNOW', out_idx==snow)
            self.__set_var('GRAUPEL', out_idx==graupel)
            self.__set_var('COD', out_idx==cod)
        
        buf = self.__memmap[out_idx==0].reshape(vars * times * lays + times, 2)
        if not (buf[:,0]==buf[:,1]).all():
            raise ValueError, "Buffer"
        
        return self.variables[key]

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
               
    def testCR(self):
        import PseudoNetCDF.testcase
        crfile = cloud_rain(PseudoNetCDF.testcase.camxfiles_paths['cloud_rain'], 4, 5)
        crfile.variables['TFLAG']
        self.assert_((crfile.variables['COD']==array([  1.25412483e+01, 1.77024829e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.38041372e+01, 1.94885385e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.41415815e+01, 1.67501605e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.34467077e+01, 1.99459922e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.25412483e+01, 1.77024829e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.38041372e+01, 1.94885385e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.41415815e+01, 1.67501605e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.34467077e+01, 1.99459922e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.25412483e+01, 1.77024829e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.38041372e+01, 1.94885385e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.41415815e+01, 1.67501605e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.34467077e+01, 1.99459922e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.96655331e+01, 2.05677104e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.14273071e+01, 2.09934115e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.21391239e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.26519203e+01, 4.96763992e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.96655331e+01, 2.05677104e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.14273071e+01, 2.09934115e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.21391239e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.26519203e+01, 4.96763992e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.96655331e+01, 2.05677104e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.14273071e+01, 2.09934115e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.21391239e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.26519203e+01, 4.96763992e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],dtype = 'f').reshape(2, 3, 4, 5)).all())
       
if __name__ == '__main__':
    unittest.main()
