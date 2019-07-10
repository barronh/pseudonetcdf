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


class ProfileTest(unittest.TestCase):
    def setUp(self):
        pass

    def testIcon(self):
        from ..cmaqfiles import icon_profile
        from ..testcase import cmaqfiles_paths
        icon_file = icon_profile(cmaqfiles_paths['icon_profile'])
        itestvals = {
            "O3": [
                3.5E-02, 3.5E-02, 4.E-02, 5.E-02, 6.E-02, 7.E-02],
            "ASO4I": [
                4.810E-03, 4.810E-03, 3.207E-03, 3.207E-03, 6.413E-04,
                3.207E-04
            ],
            "NUMATKN": [
                1.478E+09, 1.470E+09, 9.774E+08, 9.655E+08, 1.920E+08,
                9.584E+07
            ],
            "SRFACC": [
                1.491E-05, 1.373E-05, 8.694E-06, 6.815E-06, 1.190E-06,
                5.697E-07
            ]
        }
        test = self.assertEqual
        for testk, testvals in itestvals.items():
            vals = icon_file.variables[testk]
            test(True, np.allclose(vals, testvals))

    def testBcon(self):
        from ..cmaqfiles import bcon_profile
        from ..testcase import cmaqfiles_paths
        btestvals = {
            "O3": [
                3.000E-02, 3.500E-02, 4.000E-02, 5.000E-02, 6.000E-02,
                7.000E-02, 3.000E-02, 3.500E-02, 4.000E-02, 5.000E-02,
                6.000E-02, 7.000E-02, 3.500E-02, 3.500E-02, 4.000E-02,
                5.000E-02, 6.000E-02, 7.000E-02, 3.500E-02, 4.000E-02,
                4.500E-02, 5.000E-02, 6.000E-02, 7.000E-02
            ],
            "ASO4I": [
                6.413E-03, 6.413E-03, 6.413E-03, 3.207E-03, 6.413E-04,
                3.207E-04, 6.413E-03, 6.413E-03, 6.413E-03, 6.413E-03,
                6.413E-04, 3.207E-04, 4.810E-03, 4.810E-03, 3.207E-03,
                3.207E-03, 6.413E-04, 3.207E-04, 9.620E-03, 6.413E-03,
                6.413E-03, 3.207E-03, 6.413E-04, 3.207E-04
            ],
            "NUMATKN": [
                1.957E+09, 1.949E+09, 1.936E+09, 9.655E+08, 1.920E+08,
                9.584E+07, 1.957E+09, 1.949E+09, 1.936E+09, 1.924E+09,
                1.920E+08, 9.584E+07, 1.478E+09, 1.470E+09, 9.774E+08,
                9.655E+08, 1.920E+08, 9.584E+07, 2.915E+09, 1.949E+09,
                1.936E+09, 9.655E+08, 1.920E+08, 9.584E+07
            ],
            "SRFACC": [
                1.776E-05, 1.658E-05, 1.439E-05, 6.815E-06, 1.190E-06,
                5.697E-07, 1.776E-05, 1.658E-05, 1.439E-05, 1.251E-05,
                1.190E-06, 5.697E-07, 1.491E-05, 1.373E-05, 8.694E-06,
                6.815E-06, 1.190E-06, 5.697E-07, 2.346E-05, 1.658E-05,
                1.439E-05, 6.815E-06, 1.190E-06, 5.697E-07
            ],
        }
        bcon_file = bcon_profile(cmaqfiles_paths['bcon_profile'])
        test = self.assertEqual
        for testk, testvals in btestvals.items():
            vals = bcon_file.variables[testk]
            test(True, np.allclose(vals[:].T.ravel(), testvals))
