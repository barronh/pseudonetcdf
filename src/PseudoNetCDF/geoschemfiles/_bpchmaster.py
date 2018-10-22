import unittest
import os

import numpy as np

from PseudoNetCDF.pncwarn import warn
from ._bpch import bpch1
from ._newbpch import bpch2
from collections import OrderedDict


class bpch(bpch1, bpch2):
    @classmethod
    def isMine(self, *args, **kwds):
        isbpch = bpch1.isMine(*args, **kwds)
        return isbpch

    def __init__(self, bpch_path, tracerinfo=None, diaginfo=None, mode='r',
                 timeslice=slice(None), noscale=False,
                 vertgrid='GEOS-5-REDUCED', nogroup=False, reader=None):
        bpch1kwds = OrderedDict(
            tracerinfo=tracerinfo, diaginfo=diaginfo, mode=mode,
            timeslice=timeslice, noscale=noscale,
            vertgrid=vertgrid, nogroup=nogroup
        )
        bpch2kwds = OrderedDict(
            nogroup=nogroup, noscale=noscale,
            vertgrid=vertgrid
        )
        quiet = reader is None
        if reader is None or reader == 'bpch1':
            reader1, kwds1 = bpch1, bpch1kwds
            reader2, kwds2 = bpch2, bpch2kwds
        else:
            reader2, kwds2 = bpch1, bpch1kwds
            reader1, kwds1 = bpch2, bpch2kwds

        try:
            reader1.__init__(self, bpch_path, **kwds1)
        except Exception as e:
            if not quiet:
                warn('Reverting to {} : {}'.format(reader2, str(e)))
            reader2.__init__(self, bpch_path, **kwds2)


class TestMemmaps(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF.testcase import geoschemfiles_paths
        self.bpchpath = geoschemfiles_paths['bpch']

    def testNCF2BPCH(self):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', 'Not scaling variables; good for direct writing'
            )
            bpchfile = bpch(self.bpchpath, noscale=True)
        outpath = self.bpchpath + '.check'
        from PseudoNetCDF.pncgen import pncgen
        pncgen(bpchfile, outpath, inmode='r',
               outmode='w', format='bpch', verbose=0)
        orig = open(self.bpchpath, 'rb').read()
        new = open(outpath, 'rb').read()
        assert(orig == new)
        os.remove(outpath)
        from PseudoNetCDF.sci_var import reduce_dim, slice_dim
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        ALD2_check = np.array(
            [1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
             1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
             2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
             3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
             3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02,
             1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
             1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
             2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
             3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
             3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02,
             1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
             1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
             2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
             3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
             3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02])\
            .reshape(ALD2.shape)
        slided_reduced_bpchfile = slice_dim(
            reduce_dim(bpchfile, 'layer,mean'), 'time,0')
        pncgen(slided_reduced_bpchfile, outpath, inmode='r',
               outmode='w', format='bpch', verbose=0)
        ALD2_check_slided_reduced = ALD2_check[0].mean(0)[None, None]
        ALD2 = slided_reduced_bpchfile.variables['IJ-AVG-$_ALD2']
        np.testing.assert_allclose(ALD2, ALD2_check_slided_reduced * 1e-9)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', 'Not scaling variables; good for direct writing'
            )
            bpchfile = bpch(outpath)
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        np.testing.assert_allclose(ALD2, ALD2_check_slided_reduced)

    def testBPCH(self):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', 'Not scaling variables; good for direct writing'
            )
            bpchfile = bpch(self.bpchpath)
        ALD2 = bpchfile.variables['IJ-AVG-$_ALD2']
        ALD2g = bpchfile.groups['IJ-AVG-$'].variables['ALD2']
        ALD2_check = np.array(
            [1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
             1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
             2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
             3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
             3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02,
             1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
             1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
             2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
             3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
             3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02,
             1.60520077e-02, 1.82803553e-02, 2.00258084e-02, 2.01461259e-02,
             1.84865110e-02, 2.49667447e-02, 2.73083989e-02, 2.87465211e-02,
             2.89694592e-02, 2.87686456e-02, 2.87277419e-02, 3.08121163e-02,
             3.22086290e-02, 3.35262120e-02, 3.41329686e-02, 3.05218045e-02,
             3.30278911e-02, 3.58164124e-02, 3.93186994e-02, 4.15412188e-02])\
            .reshape(ALD2.shape)
        np.testing.assert_allclose(ALD2, ALD2_check)
        np.testing.assert_allclose(ALD2g, ALD2_check)
        np.testing.assert_allclose(bpchfile.variables['hyai'], np.array(
            [0.0, 0.04804826, 6.593752, 13.1348, 19.61311, 26.09201, 32.57081,
             38.98201, 45.33901, 51.69611, 58.05321, 64.36264, 70.62198,
             78.83422, 89.09992, 99.36521, 109.1817, 118.9586, 128.6959,
             142.91, 156.26, 169.609, 181.619, 193.097, 203.259, 212.15,
             218.776, 223.898, 224.363, 216.865, 201.192, 176.93, 150.393,
             127.837, 108.663, 92.36572, 78.51231, 56.38791, 40.17541,
             28.36781, 19.7916, 9.292942, 4.076571, 1.65079, 0.6167791,
             0.211349, 0.06600001, 0.01], 'f'))
        np.testing.assert_allclose(bpchfile.variables['hybi'], np.array(
            [1.0, 0.984952, 0.963406, 0.941865, 0.920387, 0.898908, 0.877429,
             0.856018, 0.8346609, 0.8133039, 0.7919469, 0.7706375, 0.7493782,
             0.721166, 0.6858999, 0.6506349, 0.6158184, 0.5810415, 0.5463042,
             0.4945902, 0.4437402, 0.3928911, 0.3433811, 0.2944031, 0.2467411,
             0.2003501, 0.1562241, 0.1136021, 0.06372006, 0.02801004,
             0.006960025, 8.175413e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'f'))
        etai_check = np.array(
            [1.01325000e+03, 9.98050662e+02, 9.82764882e+02, 9.67479511e+02,
             9.52195238e+02, 9.36910541e+02, 9.21625744e+02, 9.06342248e+02,
             8.91059167e+02, 8.75776287e+02, 8.60493406e+02, 8.45211087e+02,
             8.29929441e+02, 8.09555669e+02, 7.84087994e+02, 7.58621022e+02,
             7.33159694e+02, 7.07698900e+02, 6.82238631e+02, 6.44053520e+02,
             6.05879758e+02, 5.67705907e+02, 5.29549900e+02, 4.91400941e+02,
             4.53269420e+02, 4.15154739e+02, 3.77070069e+02, 3.39005328e+02,
             2.88927351e+02, 2.45246173e+02, 2.08244245e+02, 1.76930008e+02,
             1.50393000e+02, 1.27837000e+02, 1.08663000e+02, 9.23657200e+01,
             7.85123100e+01, 5.63879100e+01, 4.01754100e+01, 2.83678100e+01,
             1.97916000e+01, 9.29294200e+00, 4.07657100e+00, 1.65079000e+00,
             6.16779100e-01, 2.11349000e-01, 6.60000100e-02, 1.00000000e-02])
        np.testing.assert_allclose(bpchfile.variables['etai_pressure'],
                                   etai_check)

    def runTest(self):
        pass


if __name__ == '__main__':
    unittest.main()
