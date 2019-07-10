import unittest
import numpy as np


np_all_close = np.testing.assert_allclose

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

