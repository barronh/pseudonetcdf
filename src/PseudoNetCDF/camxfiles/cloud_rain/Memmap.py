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

# Distribution packages
import unittest
import struct
from warnings import warn

# Site-Packages
from numpy import zeros, array, memmap
import numpy as np

# This Package modules
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable
from PseudoNetCDF.sci_var import PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

# for use in identifying uncaught nan
listnan = struct.unpack('>f', b'\xff\xc0\x00\x00')[0]
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

    def __init__(self, rf, rows=None, cols=None):
        f = open(rf, 'rb')
        f.seek(0, 2)
        flen = f.tell()
        offset = struct.unpack('>i', open(rf, 'rb').read(4))[0] + 8
        self.__memmap = memmap(rf, dtype='>f', mode='r', offset=offset)
        cldhdrlen = offset - 20
        line1fmt = '>i%dciiii' % cldhdrlen
        ncols, nrows, nlays = struct.unpack(line1fmt,
                                            open(rf, 'rb').read(offset))[-4:-1]
        self.createDimension('COL', ncols)
        self.createDimension('ROW', nrows)
        self.createDimension('LAY', nlays)
        mydt = '>i4,S%d,>i4,>i4,>i4,>i4,>i4,>f4,>i4' % cldhdrlen
        header = np.fromfile(rf, dtype=mydt, count=1)[0]
        self.FILEDESC = ''.join(header[1].decode())
        self.STIME, self.SDATE = header.tolist()[-2:]
        if self.SDATE < 10000:
            self.SDATE += 2000000
        if (((ncols != cols and cols is not None) or
             (rows != rows and rows is not None))):
            warn(('Files says cols = %d, rows = %d, and lays = %d; ' +
                  'you said cols = %d and rows = %d') % (
                ncols, nrows, nlays, cols, rows))

        self.createDimension('DATE-TIME', 2)

        datasize = (flen - offset)
        # Try 5 first (contemporary)
        # Try 3 second (old)
        # end on 5 as a failsafe
        for nvars in [5, 3, 5]:
            timesize = (nvars * nlays * (nrows * ncols + 2) * 4 + 16)
            if (datasize % timesize) == 0:
                break
        else:
            warn(
                'File appears incomplete using 3 (v4.2) or 5 ' +
                'variables (>=v4.3); expected to fail'
            )

        if nvars < 5:
            self.VERSION = '<4.3'
            varkeys = ['CLOUD', 'PRECIP', 'COD', 'TFLAG']
        else:
            self.VERSION = '>=4.3'
            varkeys = ['CLOUD', 'RAIN', 'SNOW', 'GRAUPEL', 'COD', 'TFLAG']

        ntimes = datasize // timesize
        self.createDimension('TSTEP', ntimes)
        self.createDimension('VAR', len(varkeys) - 1)

        self.NVARS = len(self.dimensions['VAR'])
        self.NLAYS = len(self.dimensions['LAY'])
        self.NROWS = len(self.dimensions['ROW'])
        self.NCOLS = len(self.dimensions['COL'])
        self.FTYPE = 1

        self.variables = PseudoNetCDFVariables(self.__var_get, varkeys)

        self.SDATE, self.STIME = self.variables['TFLAG'][0, 0, :]

    def __set_var(self, key, vals_idx):
        times = len(self.dimensions['TSTEP'])
        lays = len(self.dimensions['LAY'])
        rows = len(self.dimensions['ROW'])
        cols = len(self.dimensions['COL'])
        vals = self.__memmap[vals_idx].reshape(times, lays, rows, cols)
        v = PseudoNetCDFVariable(self, key, 'f',
                                 ('TSTEP', 'LAY', 'ROW', 'COL'),
                                 values=vals)
        v.units = {'COD': 'None'}.get(key, 'g/m**3')
        v.long_name = key
        v.var_desc = key
        self.variables[key] = v

    def __var_get(self, key):
        times = len(self.dimensions['TSTEP'])
        rows = len(self.dimensions['ROW'])
        cols = len(self.dimensions['COL'])
        lays = len(self.dimensions['LAY'])
        vars = len(list(self.variables.keys())) - 1
        hour = 1
        date = 2
        cloud = 3
        rain = 4
        snow = 5
        graupel = 6
        cod = 7
        # stagger = 8
        out_idx = zeros(self.__memmap.shape, dtype='b')
        out_idx.reshape(times, lays * vars *
                        (rows * cols + 2) + 4)[:, 1] = hour
        out_idx.reshape(times, lays * vars *
                        (rows * cols + 2) + 4)[:, 2] = date

        dateblock = self.__memmap[out_idx == date].view('>i')
        hourblock = self.__memmap[out_idx == hour]
        nvars = len(self.dimensions['VAR'])
        self.variables['TFLAG'] = ConvertCAMxTime(dateblock, hourblock, nvars)

        newshape1 = (times, lays * vars * (rows * cols + 2) + 4)
        newshape2 = (times, lays, vars, rows * cols + 2)
        newshape3 = (times, lays, vars, rows, cols)
        val_shape = out_idx.reshape(*newshape1)[:, 4:]
        val_shape = val_shape.reshape(*newshape2)[:, :, :, 1:-1]
        val_shape = val_shape.reshape(*newshape3)
        if self.VERSION == '<4.3':
            val_shape[:, :, 0, :, :] = cloud
            val_shape[:, :, 1, :, :] = rain
            val_shape[:, :, 2, :, :] = cod
            self.__set_var('CLOUD', out_idx == cloud)
            self.__set_var('PRECIP', out_idx == rain)
            self.__set_var('COD', out_idx == cod)
        else:
            val_shape[:, :, 0, :, :] = cloud
            val_shape[:, :, 1, :, :] = rain
            val_shape[:, :, 2, :, :] = snow
            val_shape[:, :, 3, :, :] = graupel
            val_shape[:, :, 4, :, :] = cod
            self.__set_var('CLOUD', out_idx == cloud)
            self.__set_var('RAIN', out_idx == rain)
            self.__set_var('SNOW', out_idx == snow)
            self.__set_var('GRAUPEL', out_idx == graupel)
            self.__set_var('COD', out_idx == cod)

        buf = self.__memmap[out_idx == 0].reshape(
            vars * times * lays + times, 2)
        if not (buf[:, 0] == buf[:, 1]).all():
            raise ValueError("Buffer")

        return self.variables[key]


class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass

    def setUp(self):
        pass

    def testCR(self):
        import PseudoNetCDF.testcase
        crfile = cloud_rain(
            PseudoNetCDF.testcase.camxfiles_paths['cloud_rain'], 4, 5)
        crfile.variables['TFLAG']
        checkv = array([1.25412483e+01, 1.77024829e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 1.38041372e+01,
                        1.94885385e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 1.41415815e+01, 1.67501605e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        1.34467077e+01, 1.99459922e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 1.25412483e+01,
                        1.77024829e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 1.38041372e+01, 1.94885385e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        1.41415815e+01, 1.67501605e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 1.34467077e+01,
                        1.99459922e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 1.25412483e+01, 1.77024829e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        1.38041372e+01, 1.94885385e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 1.41415815e+01,
                        1.67501605e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 1.34467077e+01, 1.99459922e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        1.96655331e+01, 2.05677104e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 2.14273071e+01,
                        2.09934115e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 2.21391239e+01, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        2.26519203e+01, 4.96763992e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 1.96655331e+01,
                        2.05677104e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 2.14273071e+01, 2.09934115e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        2.21391239e+01, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 2.26519203e+01,
                        4.96763992e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 1.96655331e+01, 2.05677104e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        2.14273071e+01, 2.09934115e+00, 0.00000000e+00,
                        0.00000000e+00, 0.00000000e+00, 2.21391239e+01,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                        0.00000000e+00, 2.26519203e+01, 4.96763992e+00,
                        0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                       dtype='f').reshape(2, 3, 4, 5)
        self.assertTrue((crfile.variables['COD'] == checkv).all())

    def testNCF2CR(self):
        import PseudoNetCDF.testcase
        from PseudoNetCDF.pncgen import pncgen
        import os
        inpath = PseudoNetCDF.testcase.camxfiles_paths['cloud_rain']
        outpath = inpath + '.check'
        infile = cloud_rain(inpath, 4, 5)
        pncgen(infile, outpath, format='camxfiles.cloud_rain')
        orig = open(inpath, 'rb').read()
        new = open(outpath, 'rb').read()
        assert(orig == new)
        os.remove(outpath)


if __name__ == '__main__':
    unittest.main()
