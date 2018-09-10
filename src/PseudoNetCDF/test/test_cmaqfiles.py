import unittest
import numpy as np


np_all_close = np.testing.assert_allclose


class IOAPITest(unittest.TestCase):
    def setUp(self):
        from datetime import datetime, timedelta
        self.testncf = self._makencf()
        self.mymeta = set([
            'time', 'time_bounds', 'latitude', 'longitude', 'latitude_bounds',
            'longitude_bounds', 'lambert_conformal_conic'
        ])
        self.myvars = self.mymeta.union(['O3'])
        self.mydims = ['TIME', 'LAY', 'ROW', 'COL', 'nv', 'tnv']
        rtime = datetime.strptime(
            '1970-01-01 00:00:00+0000', '%Y-%m-%d %H:%M:%S%z')
        self.mytimes = np.array([rtime + timedelta(hours=i)
                                 for i in range(24)])
        self.testf = self._makencf()

    def _makencf(self):
        from PseudoNetCDF.cmaqfiles import ioapi_base
        minef = ioapi_base()

        minef.createDimension('TSTEP', 24)
        minef.createDimension('DATE-TIME', 2)
        minef.createDimension('LAY', 3)
        minef.createDimension('VAR', 2)
        minef.createDimension('ROW', 5)
        minef.createDimension('COL', 6)

        minef.EXEC_ID = "mcip ".ljust(80)
        minef.FTYPE = 1
        minef.CDATE = 2017069
        minef.CTIME = 144432
        minef.WDATE = 2017069
        minef.WTIME = 144432
        minef.SDATE = 2011001
        minef.STIME = 0
        minef.TSTEP = 10000
        minef.NTHIK = 1
        minef.NCOLS = 6
        minef.NROWS = 5
        minef.NLAYS = 3
        minef.NVARS = 2
        minef.GDTYP = 2
        minef.P_ALP = 33.
        minef.P_BET = 45.
        minef.P_GAM = -97.
        minef.XCENT = -97.
        minef.YCENT = 40.
        minef.XORIG = -2736000.
        minef.YORIG = -2088000.
        minef.XCELL = 36000.
        minef.YCELL = 36000.
        minef.VGTYP = 7
        minef.VGTOP = 5000.
        minef.VGLVLS = np.array([1., 0.9975, 0.995, 0.9925], dtype='f')
        minef.GDNAM = "METCRO_36US1_CRO"
        minef.UPNAM = "METCRO          "

        o3 = minef.createVariable(
            'O3', 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), units='ppbV'
        )
        no2 = minef.createVariable(
            'NO2', 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), units='ppbV'
        )

        o3[:] = np.arange(24 * 3 * 5 * 6).reshape(24, 3, 5, 6)
        no2[:] = np.arange(24 * 3 * 5 * 6).reshape(24, 3, 5, 6) * .05
        return minef

    def testTFLAG(self):
        testf = self.testf.copy()
        testf.updatetflag()
        test = self.assertEqual
        test(2011001, testf.variables['TFLAG'][0, 0, 0])
        test(0, testf.variables['TFLAG'][0, 0, 1])
        test(testf.SDATE, testf.variables['TFLAG'][0, 0, 0])
        test(testf.STIME, testf.variables['TFLAG'][0, 0, 1])
        test(testf.NVARS, testf.variables['TFLAG'].shape[1])

    def testSubsetVariables(self):
        testf = self.testf.copy()
        testf.updatetflag()
        test = self.assertEqual
        outf = testf.subsetVariables(['O3'])
        test(1, outf.variables['TFLAG'].shape[1])
        test(outf.NVARS, outf.variables['TFLAG'].shape[1])
        test('O3'.ljust(16), getattr(outf, 'VAR-LIST'))

    def _trimdim(self, dim):
        reff = self.testf.copy()
        reff.updatetflag()
        newf = reff.sliceDimensions(**{dim: slice(1, -1)})
        test = self.assertEqual
        attrkey = 'N{}S'.format(dim)
        test(getattr(reff, attrkey) - 2, getattr(newf, attrkey))
        test(len(newf.dimensions[dim]), getattr(newf, attrkey))
        return reff, newf

    def testSliceRow(self):
        reff, newf = self._trimdim('ROW')
        self.assertEqual(newf.YORIG, reff.YORIG + reff.YCELL)

    def testSliceCol(self):
        reff, newf = self._trimdim('COL')
        self.assertEqual(newf.XORIG, reff.XORIG + reff.XCELL)

    def testSliceLay(self):
        reff, newf = self._trimdim('LAY')
        self.assertEqual(True, np.allclose(newf.VGLVLS, reff.VGLVLS[1:-1]))

    def _meandim(self, dim):
        reff = self.testf.copy()
        reff.updatetflag()
        newf = reff.applyAlongDimensions(**{dim: 'mean'})
        test = self.assertEqual
        attrkey = 'N{}S'.format(dim)
        test(1, getattr(newf, attrkey))
        test(len(newf.dimensions[dim]), getattr(newf, attrkey))
        return reff, newf

    def testMeanRow(self):
        reff, newf = self._meandim('ROW')

    def testMeanCol(self):
        reff, newf = self._meandim('COL')

    def testMeanLay(self):
        reff, newf = self._meandim('LAY')
        test = self.assertEqual
        newl = reff.VGLVLS[:-1].mean()
        newu = reff.VGLVLS[1:].mean()
        test(True, np.allclose([newl, newu], newf.VGLVLS))

    def testInterpSigma(self):
        reff = self.testf.copy()
        reff.updatetflag()
        # Simple resample
        newf = reff.interpSigma(
            reff.VGLVLS[1:], vgtop=reff.VGTOP, interptype='conserve'
        )
        test = np.allclose(
            reff.variables['O3'][:, 1:], newf.variables['O3'][:]
        )
        self.assertEqual(True, test)
        newf = reff.interpSigma(
            np.convolve([.5, .5], reff.VGLVLS, mode='valid'),
            vgtop=reff.VGTOP, interptype='conserve'
        )
        test = np.allclose(
            0.5 * reff.variables['O3'][:, 1:] +
            0.5 * reff.variables['O3'][:, :-1],
            newf.variables['O3'][:]
        )
