__all__ = ['landuse']
__doc__ =  """
.. _Memmap
:mod:`Memmap` -- landuse Memmap interface
 ============================================ 

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              landuse files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL = "$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate =  "$LastChangedDate$"
RevisionNum = "$LastChangedRevision$"
ChangedBy =   "$LastChangedBy: svnbarronh $"
__version__ =  RevisionNum

#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros, array, where, memmap, newaxis, dtype, nan

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff, timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile, Int2Asc
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan = struct.unpack('>f', '\xff\xc0\x00\x00')[0]
checkarray = zeros((1, ), 'f')
checkarray[0] = listnan
array_nan = checkarray[0]
    
class landuse(PseudoNetCDFFile):
    """
    landuse provides a PseudoNetCDF interface for CAMx
    landuse files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> landuse_path =  'camx_landuse.bin'
        >>> rows, cols =  65, 83
        >>> landusefile =  landuse(landuse_path, rows, cols)
        >>> landusefile.variables.keys()
        ['TFLAG', 'FLAND', 'TOPO']
        >>> tflag =  landusefile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0, 0, :]
        array([2005185,      0])
        >>> tflag[-1, 0, :]
        array([2005185, 240000])
        >>> v =  landusefile.variables['FLAND']
        >>> v.dimensions
        ('LANDUSE', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> landusefile.dimensions
        {'LANDUSE': 11, 'ROW': 65, 'COL': 83, 'VAR': 2}
    """
    
    def __init__(self, rf, rows, cols, mode = 'r'):
        self.__mode = mode
        self._rffile = OpenRecordFile(rf)
        self._rffile.infile.seek(0, 2)
        self.__rf = rf
        rflen = self._rffile.infile.tell()
        self._rffile._newrecord(0)
        
        self.createDimension('ROW', rows)
        self.createDimension('COL', cols)
        first_line, =  self._rffile.read('8s')
        if first_line ==  'LUCAT11 ':
            self.createDimension('LANDUSE', 11)
            self._newstyle =  True
        elif first_line ==  'LUCAT26 ':
            self.createDimension('LANDUSE', 26)
            self._newstyle =  True
        else:
            self.createDimension('LANDUSE', 11)
            self._newstyle =  False
        nland = len(self.dimensions['LANDUSE'])
        nrows = len(self.dimensions['ROW'])
        ncols = len(self.dimensions['COL'])
        if self._newstyle:
            self.__fland_dtype = dtype(dict(names = ['SPAD1', 'KEY', 'EPAD1', 'SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '8>S', '>i', '>i', '(%d, %d, %d)>f' % (nland, nrows, ncols), '>i']))
            self.__other_dtype = dtype(dict(names = ['SPAD1', 'KEY', 'EPAD1', 'SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '8>S', '>i', '>i', '(%d, %d)>f' % (nrows, ncols), '>i']))
        else:
            self.__fland_dtype = dtype(dict(names = ['SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '(%d, %d, %d)>f' % (nland, nrows, ncols), '>i']))
            self.__other_dtype = dtype(dict(names = ['SPAD2', 'DATA', 'EPAD2'], formats = ['>i', '(%d, %d)>f' % (nrows, ncols), '>i']))
            
        self.__addvars()
        if self._newstyle:
            self.__keys = [first_line]
            
        else:
            self.__keys = ['LUCAT11']

    def __addvars(self):
        nrows =  len(self.dimensions['ROW'])
        ncols =  len(self.dimensions['COL'])
        nland =  len(self.dimensions['LANDUSE'])
        self._rffile.infile.seek(0, 2)
        rflen = self._rffile.infile.tell()
        fland_dtype = self.__fland_dtype
        other_dtype = self.__other_dtype
        nfland = fland_dtype.itemsize
        nfland1opt = nfland + other_dtype.itemsize
        nfland2opt = nfland + other_dtype.itemsize * 2
        if rflen == nfland:
            file_dtype = dtype(dict(names = ['FLAND'], formats = [fland_dtype]))
        elif rflen == nfland1opt:
            file_dtype = dtype(dict(names = ['FLAND', 'VAR1'], formats = [fland_dtype, other_dtype]))
        elif rflen == nfland2opt:
            file_dtype = dtype(dict(names = ['FLAND', 'LAI', 'TOPO'], formats = [fland_dtype, other_dtype, other_dtype]))
        else:
            raise IOError('File size is expected to be %d or %d; was %d' % (nfland, nfland1opt, nfland2opt))
        
        data = memmap(self.__rf, mode = self.__mode, dtype = file_dtype, offset = 0)
        if not self._newstyle:
            varkeys = ['FLAND', 'TOPO']
        else:
            varkeys = [data[k]['KEY'][0].strip() for k in file_dtype.names]
        
        for varkey, dkey in zip(varkeys, file_dtype.names):
            var = self.createVariable(varkey, 'f', {'FLAND': ('LANDUSE', 'ROW', 'COL')}.get(dkey, ('ROW', 'COL')))
            var[:] = data[dkey]['DATA']
            var.var_desc = varkey.ljust(16)
            var.units = {'FLAND': 'Fraction'}.get(dkey, '')
            var.long_name = varkey.ljust(16)
        
        

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
               
    def testLU(self):
        from numpy import testing
        import PseudoNetCDF.testcase
        aassert = testing.assert_array_almost_equal
        lufile = landuse(PseudoNetCDF.testcase.camxfiles_paths['landuse'], 4, 5)
        
        aassert(lufile.variables['LUCAT11'], array([  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 6.25000000e-02, 0.00000000e+00, 6.25000000e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 6.25000000e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 9.37500000e-01, 1.00000000e+00, 9.37500000e-01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 9.37500000e-01, 6.87500000e-01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 6.25000000e-02, 2.50000000e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00], dtype = 'f').reshape(11, 4, 5))

if __name__ ==  '__main__':
    unittest.main()
