import unittest
import numpy as np


np_all_close = np.testing.assert_allclose


class testsonde(unittest.TestCase):
    def testSonde(self):
        from ..woudcfiles import woudcsonde
        from ..testcase import woudcfiles_paths
        wf = woudcsonde(woudcfiles_paths['sonde'])
        pv = wf.variables['Pressure']
        np_all_close(
            pv[0], [
                1013.39, 950.48, 850.84, 750.23, 650.28, 550.28, 450.17,
                350.18, 250.15, 150.17, 50., 7.46
            ]
        )
        assert (pv.units.strip() == 'hPa')
        assert (pv.long_name.strip() == 'Pressure')
        wqf = wf.avgSigma(
            np.array([1, 0]), 0
        )
        newo3 = wqf.variables['O3']
        np_all_close(newo3, [[0.0596650205552578]])
        assert (newo3.units.strip() == 'ppm')
