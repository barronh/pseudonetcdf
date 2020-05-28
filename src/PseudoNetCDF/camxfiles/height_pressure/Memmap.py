__all__ = ['height_pressure']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- height/pressure Memmap interface
=================================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              height/pressure files.  See PseudoNetCDFFile
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

# Distribution packages
import unittest
import struct

# Site-Packages
from numpy import zeros, array, memmap

# This Package modules
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable
from PseudoNetCDF.sci_var import PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

# for use in identifying uncaught nan
listnan = struct.unpack('>f', b'\xff\xc0\x00\x00')[0]
checkarray = zeros((1,), 'f')
checkarray[0] = listnan
array_nan = checkarray[0]


class height_pressure(PseudoNetCDFFile):
    """
    height_pressure provides a PseudoNetCDF interface for CAMx
    height_pressure files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).

    ex:
        >>> height_pressure_path = 'camx_height_pressure.bin'
        >>> rows,cols = 65,83
        >>> hpf = height_pressure(height_pressure_path,rows,cols)
        >>> hpf.variables.keys()
        ['TFLAG', 'HGHT', 'PRES']
        >>> v = hpf.variables['V']
        >>> tflag = hpf.variables['TFLAG']
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
        >>> hpf.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """

    id_fmt = 'fi'
    data_fmt = 'f'

    def __init__(self, rf, rows=None, cols=None):
        self.__memmap = memmap(rf, '>f', 'r', offset=0)
        rowsXcols = self.__memmap[0].view('i') // 4 - 2
        record_length = rowsXcols + 4
        records = self.__memmap.size // record_length

        times = self.__memmap.reshape(records, record_length)[:, 1:3]
        self.STIME, self.SDATE = times[0]
        for i, (t, d) in enumerate(times):
            if (t, d) != (self.STIME, self.SDATE):
                break
        self.SDATE = self.SDATE.view('i')
        self.createDimension('LAY', i / 2)
        self.createDimension('TSTEP', times.shape[0] / i)

        if rows is None and cols is None:
            rows = rowsXcols
            cols = 1
        elif rows is None:
            rows = rowsXcols / cols
        elif cols is None:
            cols = rowsXcols / rows
        else:
            if cols * rows != rowsXcols:
                raise ValueError(("The product of cols (%d) and rows (%d) " +
                                  "must equal cells (%d)") % (
                                 cols, rows, rowsXcols))

        self.createDimension('ROW', rows)
        self.createDimension('COL', cols)
        self.createDimension('DATE-TIME', 2)
        self.createDimension('VAR', 2)

        self.NROWS = rows
        self.NCOLS = cols
        self.NLAYS = len(self.dimensions['LAY'])
        self.NVARS = 2
        self.NTHIK = 1

        setattr(self, 'VAR-LIST', 'HGHT'.ljust(16) + 'PRES'.ljust(16))

        self.variables = PseudoNetCDFVariables(
            self.__var_get, ['HGHT', 'PRES', 'TFLAG'])

    def __var_get(self, key):
        lays = len(self.dimensions['LAY'])
        times = len(self.dimensions['TSTEP'])
        rows = len(self.dimensions['ROW'])
        cols = len(self.dimensions['COL'])
        hght = 1
        pres = 2
        time = 3
        date = 4
        out_idx = zeros(self.__memmap.shape, dtype='b').reshape(
            times, lays, 2, rows * cols + 4)
        out_idx[:, :, 0, 3:-1] = hght
        out_idx[:, :, 1, 3:-1] = pres
        out_idx[:, :, :, 1] = time
        out_idx[:, :, :, 2] = date
        out_idx = out_idx.ravel()
        buf = self.__memmap[out_idx == 0].reshape(lays * 2 * times, 2)
        if not (buf[:, 0] == buf[:, 1]).all():
            raise ValueError("Buffer")
        tmpvals = self.__memmap[out_idx == 1].reshape(times, lays, rows, cols)
        v = self.variables['HGHT'] = PseudoNetCDFVariable(
            self, 'HGHT', 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), values=tmpvals)
        v.units = 'm'
        v.long_name = 'HGHT'.ljust(16)
        v.var_desc = 'Top Height'
        tmpvals = self.__memmap[out_idx == 2].reshape(times, lays, rows, cols)
        v = self.variables['PRES'] = PseudoNetCDFVariable(
            self, 'PRES', 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), values=tmpvals)
        v.units = 'hPA'
        v.long_name = 'PRES'.ljust(16)
        v.var_desc = 'Pressure at center'
        tslice = slice(None, None, len(self.dimensions['LAY']) * 2)
        dateblock = self.__memmap[out_idx == 4][tslice].view('>i')
        hourblock = self.__memmap[out_idx == 3][tslice]
        self.variables['TFLAG'] = ConvertCAMxTime(dateblock, hourblock,
                                                  len(self.dimensions['VAR']))

        return self.variables[key]


class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass

    def setUp(self):
        pass

    def testHP(self):
        import PseudoNetCDF.testcase
        hpfile = height_pressure(
            PseudoNetCDF.testcase.camxfiles_paths['height_pressure'], 4, 5)
        hpfile.variables['TFLAG']
        checkv = array([3.38721924e+01, 3.40657959e+01, 3.41392822e+01,
                        3.42358398e+01, 3.42543945e+01, 3.38868408e+01,
                        3.40622559e+01, 3.42358398e+01, 3.44768066e+01,
                        3.46112061e+01, 3.37558594e+01, 3.39323730e+01,
                        3.42663574e+01, 3.46854248e+01, 3.48144531e+01,
                        3.39472656e+01, 3.41900635e+01, 3.46160889e+01,
                        3.48209229e+01, 3.47874756e+01, 6.78652344e+01,
                        6.82532959e+01, 6.84020996e+01, 6.85950928e+01,
                        6.86304932e+01, 6.78945312e+01, 6.82465820e+01,
                        6.85941162e+01, 6.90783691e+01, 6.93474121e+01,
                        6.76313477e+01, 6.79859619e+01, 6.86558838e+01,
                        6.94960938e+01, 6.97552490e+01, 6.80159912e+01,
                        6.85028076e+01, 6.93570557e+01, 6.97674561e+01,
                        6.97009277e+01, 1.01980713e+02, 1.02563843e+02,
                        1.02787109e+02, 1.03077759e+02, 1.03131104e+02,
                        1.02022949e+02, 1.02553101e+02, 1.03076904e+02,
                        1.03804565e+02, 1.04208984e+02, 1.01628662e+02,
                        1.02162842e+02, 1.03169922e+02, 1.04433838e+02,
                        1.04823120e+02, 1.02206909e+02, 1.02940186e+02,
                        1.04224609e+02, 1.04841797e+02, 1.04740723e+02,
                        3.38721924e+01, 3.40657959e+01, 3.41392822e+01,
                        3.42358398e+01, 3.42543945e+01, 3.38868408e+01,
                        3.40622559e+01, 3.42358398e+01, 3.44768066e+01,
                        3.46112061e+01, 3.37558594e+01, 3.39323730e+01,
                        3.42663574e+01, 3.46854248e+01, 3.48144531e+01,
                        3.39472656e+01, 3.41900635e+01, 3.46160889e+01,
                        3.48209229e+01, 3.47874756e+01, 6.78652344e+01,
                        6.82532959e+01, 6.84020996e+01, 6.85950928e+01,
                        6.86304932e+01, 6.78945312e+01, 6.82465820e+01,
                        6.85941162e+01, 6.90783691e+01, 6.93474121e+01,
                        6.76313477e+01, 6.79859619e+01, 6.86558838e+01,
                        6.94960938e+01, 6.97552490e+01, 6.80159912e+01,
                        6.85028076e+01, 6.93570557e+01, 6.97674561e+01,
                        6.97009277e+01, 1.01980713e+02, 1.02563843e+02,
                        1.02787109e+02, 1.03077759e+02, 1.03131104e+02,
                        1.02022949e+02, 1.02553101e+02, 1.03076904e+02,
                        1.03804565e+02, 1.04208984e+02, 1.01628662e+02,
                        1.02162842e+02, 1.03169922e+02, 1.04433838e+02,
                        1.04823120e+02, 1.02206909e+02, 1.02940186e+02,
                        1.04224609e+02, 1.04841797e+02, 1.04740723e+02],
                       dtype='f').reshape(2, 3, 4, 5)
        self.assertTrue((hpfile.variables['HGHT'] == checkv).all())

    def testNCF2HP(self):
        import PseudoNetCDF.testcase
        from PseudoNetCDF.pncgen import pncgen
        import os
        inpath = PseudoNetCDF.testcase.camxfiles_paths['height_pressure']
        outpath = inpath + '.check'
        infile = height_pressure(inpath, 4, 5)
        pncgen(infile, outpath, format='camxfiles.height_pressure')
        orig = open(inpath, 'rb').read()
        new = open(outpath, 'rb').read()
        assert(orig == new)
        os.remove(outpath)


if __name__ == '__main__':
    unittest.main()
