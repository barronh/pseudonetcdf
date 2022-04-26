import unittest
import numpy as np
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariables
from PseudoNetCDF import PseudoNetCDFVariable, pncopen
from . import requires_basemap, requires_pyproj
from PseudoNetCDF.pncwarn import warn


np_all_close = np.testing.assert_allclose


class PseudoNetCDFFileTest(unittest.TestCase):
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

    def _makencf(self):
        from numpy import arange
        tncf = PseudoNetCDFFile()

        tncf.createDimension('TSTEP', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)
        tncf.createDimension('nv', 4)
        tncf.createDimension('tnv', 2)
        tncf.str_one = '1'
        tncf.int_two = 2
        tncf.float_threeptfive = 3.5
        tncf.Conventions = 'CF-1.6'
        o3 = tncf.createVariable('O3', 'f', ('TSTEP', 'LAY', 'ROW', 'COL'))

        o3[:] = arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)
        o3.units = 'ppbv'
        o3.grid_mapping = 'lambert_conformal_conic'
        time = tncf.createVariable('time', 'd', ('TSTEP',))
        time.long_name = 'time'
        time.units = 'hours since 1970-01-01 00:00:00+0000'
        time[:] = np.arange(24)

        timeb = tncf.createVariable('time_bounds', 'd', ('TSTEP', 'tnv'))
        timeb.long_name = 'time_bounds'
        timeb.units = 'hours since 1970-01-01 00:00:00+0000'
        timeb[:, 0] = np.arange(0, 24)
        timeb[:, 1] = np.arange(1, 25)

        crs = tncf.createVariable('lambert_conformal_conic', 'i', ())
        crs.grid_mapping_name = 'lambert_conformal_conic'
        crs.standard_parallel = np.array([30., 45.])
        crs.longitude_of_central_meridian = -97.
        crs.latitude_of_projection_origin = 40.
        crs.false_northing = 1620000.
        crs.false_easting = 2412000.
        crs.semi_major_axis = 6371000.
        crs.semi_minor_axis = 6371000.
        lat = tncf.createVariable('latitude', 'f', ('ROW', 'COL'))
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lon = tncf.createVariable('longitude', 'f', ('ROW', 'COL'))
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'
        latb = tncf.createVariable(
            'latitude_bounds', 'f', ('ROW', 'COL', 'nv'))
        latb.long_name = 'latitude_bounds'
        latb.units = 'degrees_north'
        lonb = tncf.createVariable(
            'longitude_bounds', 'f', ('ROW', 'COL', 'nv'))
        lonb.long_name = 'longitude_bounds'
        lonb.units = 'degrees_east'
        lon[:] = [
            [-120.21161038333193, -120.21160114763147, -120.21159191193058,
             -120.21158267622918, -120.21157344052737, -120.21156420482505],
            [-120.21161271536134, -120.21160347966001, -120.21159424395826,
             -120.21158500825604, -120.21157577255335, -120.21156653685021],
            [-120.21161504739118, -120.21160581168901, -120.21159657598642,
             -120.21158734028334, -120.2115781045798, -120.2115688688758],
            [-120.21161737942151, -120.21160814371851, -120.21159890801503,
             -120.21158967231109, -120.21158043660672, -120.21157120090189],
            [-120.21161971145229, -120.21161047574842, -120.21160124004409,
             -120.21159200433934, -120.21158276863409, -120.21157353292838]]

        lat[:] = [
            [22.748507533242535, 22.748509683865187, 22.74851183448702,
             22.74851398510802, 22.748516135728206, 22.748518286347593],
            [22.74851605050742, 22.748518201130356, 22.748520351752475,
             22.74852250237377, 22.74852465299425, 22.748526803613903],
            [22.74852456777256, 22.748526718395773, 22.748528869018187,
             22.748531019639763, 22.748533170260536, 22.74853532088048],
            [22.748533085037966, 22.74853523566145, 22.74853738628417,
             22.748539536906023, 22.7485416875271, 22.748543838147327],
            [22.74854160230359, 22.74854375292739, 22.748545903550376,
             22.748548054172538, 22.748550204793883, 22.7485523554144]]
        lonb[:] = [
            [[-120.21161038333193, -120.21160114763147,
              -120.21160347966001, -120.21161271536134],
             [-120.21160114763147, -120.21159191193058,
              -120.21159424395826, -120.21160347966001],
             [-120.21159191193058, -120.21158267622918,
              -120.21158500825604, -120.21159424395826],
             [-120.21158267622918, -120.21157344052737,
              -120.21157577255335, -120.21158500825604],
             [-120.21157344052737, -120.21156420482505,
              -120.21156653685021, -120.21157577255335],
             [-120.21156420482505, -120.2115549691223,
              -120.2115573011466, -120.21156653685021]],
            [[-120.21161271536134, -120.21160347966001,
              -120.21160581168901, -120.21161504739118],
             [-120.21160347966001, -120.21159424395826,
              -120.21159657598642, -120.21160581168901],
             [-120.21159424395826, -120.21158500825604,
              -120.21158734028334, -120.21159657598642],
             [-120.21158500825604, -120.21157577255335,
              -120.2115781045798, -120.21158734028334],
             [-120.21157577255335, -120.21156653685021,
              -120.2115688688758, -120.2115781045798],
             [-120.21156653685021, -120.2115573011466,
              -120.21155963317135, -120.2115688688758]],
            [[-120.21161504739118, -120.21160581168901,
              -120.21160814371851, -120.21161737942151],
             [-120.21160581168901, -120.21159657598642,
              -120.21159890801503, -120.21160814371851],
             [-120.21159657598642, -120.21158734028334,
              -120.21158967231109, -120.21159890801503],
             [-120.21158734028334, -120.2115781045798,
              -120.21158043660672, -120.21158967231109],
             [-120.2115781045798, -120.2115688688758,
              -120.21157120090189, -120.21158043660672],
             [-120.2115688688758, -120.21155963317135,
              -120.21156196519657, -120.21157120090189]],
            [[-120.21161737942151, -120.21160814371851,
              -120.21161047574842, -120.21161971145229],
             [-120.21160814371851, -120.21159890801503,
              -120.21160124004409, -120.21161047574842],
             [-120.21159890801503, -120.21158967231109,
              -120.21159200433934, -120.21160124004409],
             [-120.21158967231109, -120.21158043660672,
              -120.21158276863409, -120.21159200433934],
             [-120.21158043660672, -120.21157120090189,
              -120.21157353292838, -120.21158276863409],
             [-120.21157120090189, -120.21156196519657,
              -120.21156429722222, -120.21157353292838]],
            [[-120.21161971145229, -120.21161047574842,
              -120.21161280777879, -120.2116220434835],
             [-120.21161047574842, -120.21160124004409,
              -120.21160357207363, -120.21161280777879],
             [-120.21160124004409, -120.21159200433934,
              -120.21159433636801, -120.21160357207363],
             [-120.21159200433934, -120.21158276863409,
              -120.21158510066192, -120.21159433636801],
             [-120.21158276863409, -120.21157353292838,
              -120.21157586495535, -120.21158510066192],
             [-120.21157353292838, -120.21156429722222,
              -120.21156662924835, -120.21157586495535]]]
        latb[:] = [
            [[22.748507533242535, 22.748509683865187,
              22.748518201130356, 22.74851605050742],
             [22.748509683865187, 22.74851183448702,
              22.748520351752475, 22.748518201130356],
             [22.74851183448702, 22.74851398510802,
              22.74852250237377, 22.748520351752475],
             [22.74851398510802, 22.748516135728206,
              22.74852465299425, 22.74852250237377],
             [22.748516135728206, 22.748518286347593,
              22.748526803613903, 22.74852465299425],
             [22.748518286347593, 22.748520436966125,
              22.748528954232754, 22.748526803613903]],
            [[22.74851605050742, 22.748518201130356,
              22.748526718395773, 22.74852456777256],
             [22.748518201130356, 22.748520351752475,
              22.748528869018187, 22.748526718395773],
             [22.748520351752475, 22.74852250237377,
              22.748531019639763, 22.748528869018187],
             [22.74852250237377, 22.74852465299425,
              22.748533170260536, 22.748531019639763],
             [22.74852465299425, 22.748526803613903,
              22.74853532088048, 22.748533170260536],
             [22.748526803613903, 22.748528954232754,
              22.748537471499613, 22.74853532088048]],
            [[22.74852456777256, 22.748526718395773,
              22.74853523566145, 22.748533085037966],
             [22.748526718395773, 22.748528869018187,
              22.74853738628417, 22.74853523566145],
             [22.748528869018187, 22.748531019639763,
              22.748539536906023, 22.74853738628417],
             [22.748531019639763, 22.748533170260536,
              22.7485416875271, 22.748539536906023],
             [22.748533170260536, 22.74853532088048,
              22.748543838147327, 22.7485416875271],
             [22.74853532088048, 22.748537471499613,
              22.748545988766764, 22.748543838147327]],
            [[22.748533085037966, 22.74853523566145,
              22.74854375292739, 22.74854160230359],
             [22.74853523566145, 22.74853738628417,
              22.748545903550376, 22.74854375292739],
             [22.74853738628417, 22.748539536906023,
              22.748548054172538, 22.748545903550376],
             [22.748539536906023, 22.7485416875271,
              22.748550204793883, 22.748548054172538],
             [22.7485416875271, 22.748543838147327,
              22.7485523554144, 22.748550204793883],
             [22.748543838147327, 22.748545988766764,
              22.748554506034104, 22.7485523554144]],
            [[22.74854160230359, 22.74854375292739,
              22.74855227019359, 22.748550119569494],
             [22.74854375292739, 22.748545903550376,
              22.748554420816852, 22.74855227019359],
             [22.748545903550376, 22.748548054172538,
              22.74855657143929, 22.748554420816852],
             [22.748548054172538, 22.748550204793883,
              22.74855872206093, 22.74855657143929],
             [22.748550204793883, 22.7485523554144,
              22.748560872681754, 22.74855872206093],
             [22.7485523554144, 22.748554506034104,
              22.748563023301763, 22.748560872681754]]]
        return tncf

    def testVal2idx(self):
        ncf = PseudoNetCDFFile()
        coorde = np.arange(9, dtype='f')
        coordc = (coorde[:-1] + coorde[1:]) / 2.
        ncf.createDimension('coord', coordc.size)
        ncf.createDimension('nv', 2)
        ncf.createVariable(
            'coord', 'f', ('coord',), values=coordc
        )
        bncf = ncf.copy()
        bncf.createVariable(
            'coord_bounds', 'f', ('coord', 'nv'),
            values=coorde.repeat(2, 0)[1:-1].reshape(-1, 2)
        )
        bncf.variables['coord'].bounds = 'coord_bounds'
        cvals = [-1, .25, 4.5, 4.99, 5, 6.75, 10]
        expectedb = np.array([0, 0, 4, 4, 5, 6, 7])
        expectedn = np.array([0, 0, 4, 4, 4, 6, 7])

        def checkvals(ncf, method, clean, compare):
            withbnds = str('coord_bounds' in ncf.variables)
            prefix = withbnds + '&' + method + '&' + clean
            idx = ncf.val2idx(
                'coord', cvals,
                method=method, clean=clean, bounds='warn'
            )
            warn(prefix + ' got: ' + repr(idx))
            warn(prefix + ' chk: ' + repr(compare))
            assert(np.ma.allclose(compare, idx))

        mw = np.ma.masked_where

        nn_mask = [0] * 7
        bn_mask = [0] * 7
        nm_mask = [1, 1, 0, 0, 0, 0, 1]
        bm_mask = [1, 0, 0, 0, 0, 0, 1]
        em_mask = [1, 1, 0, 1, 1, 1, 1]
        checkvals(ncf, 'nearest', 'none', mw(nn_mask, expectedn))
        checkvals(ncf, 'bounds', 'none', mw(bn_mask, expectedb))
        checkvals(ncf, 'nearest', 'mask', mw(nm_mask, expectedn))
        checkvals(ncf, 'bounds', 'mask', mw(bm_mask, expectedb))

        checkvals(bncf, 'nearest', 'none', mw(nn_mask, expectedn))
        checkvals(bncf, 'bounds', 'none', mw(bn_mask, expectedb))
        checkvals(bncf, 'nearest', 'mask', mw(nm_mask, expectedn))
        checkvals(bncf, 'bounds', 'mask', mw(bm_mask, expectedb))

        checkvals(bncf, 'exact', 'mask', mw(em_mask, expectedb))

    def testCopyVariable(self):
        tncf = self.testncf
        var = tncf.copyVariable(
            tncf.variables['O3'], key='O3_PPB', withdata=True)
        self.assertEqual(True, (var[:] == tncf.variables['O3_PPB']).all())

    def testSubsetVariables(self):
        tncf = self.testncf.copy()
        var = tncf.copyVariable(
            tncf.variables['O3'], key='O3_PPB', withdata=True)
        var[:] *= 1e3
        var = tncf.copyVariable(
            tncf.variables['O3'], key='O3_PPT', withdata=True)
        var[:] *= 1e6
        sncf = tncf.subsetVariables(['O3_PPT'])
        self.assertEqual(len(sncf.variables), 1)
        self.assertEqual(set(sncf.variables), set(['O3_PPT']))
        sncf = tncf.subsetVariables(['O3_PPT'], exclude=True)
        self.assertEqual(len(sncf.variables), len(
            self.myvars.union(['O3_PPB'])))
        self.assertEqual(set(sncf.variables), self.myvars.union(['O3_PPB']))

    def testRenameVariables(self):
        tncf = self.testncf
        sncf = tncf.renameVariables(O3='O3_PPM')
        self.assertEqual(len(sncf.variables), len(self.myvars))
        self.assertEqual(set(sncf.variables), self.mymeta.union(['O3_PPM']))

    def testRenameDimensions(self):
        tncf = self.testncf
        sncf = tncf.renameDimensions(TSTEP='TIME')
        self.assertEqual(len(sncf.dimensions), len(tncf.dimensions))
        self.assertEqual(set(sncf.dimensions), set(self.mydims))

    def testSliceDimensionInt(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(TSTEP=0)
        self.assertEqual(len(sncf.dimensions['TSTEP']), 1)
        self.assertEqual(
            True, (sncf.variables['O3'][:] == o3[0]).all())

    def testSliceDimensionList(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(TSTEP=[0])
        self.assertEqual(len(sncf.dimensions['TSTEP']), 1)
        self.assertEqual(
            True, (sncf.variables['O3'][:] == o3[0]).all())
        sncf = tncf.sliceDimensions(TSTEP=[0, 8])
        self.assertEqual(len(sncf.dimensions['TSTEP']), 2)
        self.assertEqual(
            True, (sncf.variables['O3'][:] == o3[[0, 8]]).all())

    def testSliceDimensionComboListInt(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(TSTEP=[0, 8], ROW=2, COL=3)
        self.assertEqual(len(sncf.dimensions['TSTEP']), 2)
        self.assertEqual(len(sncf.dimensions['ROW']), 1)
        self.assertEqual(len(sncf.dimensions['COL']), 1)
        self.assertEqual(True, (sncf.variables['O3'][:] == o3[[
                         0, 8], :, 2, 3][:, :, None, None]).all())

    def testSliceDimensionMultiArray(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        i = np.arange(4)
        sncf = tncf.sliceDimensions(TSTEP=i, LAY=i, ROW=i, COL=i)
        self.assertEqual(len(sncf.dimensions['POINTS']), 4)
        self.assertEqual(
            True, (sncf.variables['O3'][:] == o3[i, i, i, i]).all())

    def testSliceDimensionSlice(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(ROW=slice(0, 3))
        self.assertEqual(len(sncf.dimensions['ROW']), 3)
        self.assertEqual(
            True, (sncf.variables['O3'][:] == o3[:, :, 0:3, :]).all())

    def testSliceDimensionMultiSlice(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        sncf = tncf.sliceDimensions(
            LAY=slice(0, 3), ROW=slice(0, 3), COL=slice(0, 3))
        self.assertEqual(len(sncf.dimensions['LAY']), 3)
        self.assertEqual(len(sncf.dimensions['ROW']), 3)
        self.assertEqual(len(sncf.dimensions['COL']), 3)
        self.assertEqual(
            True, (sncf.variables['O3'][:] == o3[:, 0:3, 0:3, 0:3]).all())

    def testApplyAlongDimensionsNamed(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        ancf = tncf.applyAlongDimensions(LAY='min')
        self.assertEqual(
            True, (ancf.variables['O3'][:] == o3.min(1, keepdims=True)).all())

    def testApplyAlongDimensionsMultiNamed(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        ancf = tncf.applyAlongDimensions(LAY='min', ROW='max')
        self.assertEqual(True, (ancf.variables['O3'][:] == o3.min(
            1, keepdims=True).max(2, keepdims=True)).all())

    def testApplyAlongDimensionsConvolve(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        ancf = tncf.applyAlongDimensions(TSTEP=lambda x: np.convolve(
            x, np.ones(2, dtype='f') / 2., mode='valid'))
        co3 = (o3[1:] + o3[:-1]) / 2
        self.assertEqual(True, (ancf.variables['O3'][:] == co3).all())

    def testApplyAlongDimensionsConvolveWithMax(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:]
        co3 = (o3[1:] + o3[:-1]) / 2
        # Testing convolution; useful for mda8

        def convf(x):
            return np.convolve(x, np.ones(2, dtype='f') / 2., mode='valid')

        ancf = tncf.applyAlongDimensions(TSTEP=convf)\
                   .applyAlongDimensions(TSTEP=np.max)
        mco3 = co3.max(0, keepdims=True)
        self.assertEqual(True, (ancf.variables['O3'][:] == mco3).all())

    def testAdd(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        tmpf = tncf + tncf
        np_all_close(tmpf.variables['O3'][:], 2 * o3)

    def testSub(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        tmpf = tncf - tncf
        np_all_close(tmpf.variables['O3'][:], o3 - o3)

    def testDiv(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        np.seterr(invalid='ignore')
        tmpf = tncf / tncf
        np_all_close(tmpf.variables['O3'][:], o3 / o3)
        np.seterr(invalid='warn')

    def testMul(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        tmpf = tncf * tncf
        np_all_close(tmpf.variables['O3'][:], o3**2)

    def testPow(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        np.seterr(over='ignore', invalid='ignore')
        tmpf = tncf ** tncf
        np_all_close(tmpf.variables['O3'][:], o3**o3)
        np.seterr(over='warn', invalid='warn')

    def testMod(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        np.seterr(divide='ignore')
        tmpf = tncf % tncf
        np_all_close(tmpf.variables['O3'][:], o3 % o3)
        np.seterr(divide='warn')

    def testXAdd(self):
        tncf = self.testncf
        o3 = tncf.variables['O3'][:].copy()
        tmpf = tncf.copy()
        print(tmpf._operator_exclude_vars)
        tmpf.setCoords(['O3'])
        print(tmpf._operator_exclude_vars)
        tmpf = tmpf + tncf
        np_all_close(tmpf.variables['O3'][:], o3)

    @requires_basemap
    def testGetMap(self):
        tncf = self.testncf
        tncf.getMap(maptype='basemap_auto')

    @requires_pyproj
    def testGetproj(self):
        tncf = self.testncf
        tncf.getproj(withgrid=False, projformat='pyproj')

    @requires_pyproj
    def testLl2xy(self):
        tncf = self.testncf
        crs = tncf.variables['lambert_conformal_conic']
        y0 = crs.false_northing
        x0 = crs.false_easting
        lon0, lat0 = tncf.xy2ll(x0, y0)
        x0t, y0t = tncf.ll2xy(lon0, lat0)
        self.assertEqual((x0, y0), (x0t, y0t))

    @requires_pyproj
    def testLl2ij(self):
        tncf = self.testncf
        y0 = 0
        x0 = 0
        lon0, lat0 = tncf.xy2ll(x0, y0)
        i0, j0 = tncf.ll2ij(lon0, lat0)
        self.assertEqual(0, i0)
        self.assertEqual(0, j0)

    def testTime2t(self):
        from datetime import datetime
        time = np.array([
            datetime.strptime('1970-01-{}+0000'.format(t), '%Y-%m-%d %H:%M%z')
            for t in [
                '01 00:20', '01 01:20', '01 01:40', '01 06:00', '02 01:00'
            ]
        ])
        tncf = self.testncf
        t = tncf.time2t(time[:])
        np_all_close(t, [0, 1, 2, 6, 23])
        t = tncf.time2t(time[:-1], ttype='bounds')
        np_all_close(t, [0, 1, 1, 6])
        t = tncf.time2t(time[:], ttype='bounds_close')
        self.assertEqual(True, np.allclose(t, [0, 1, 1, 6, 23]))
        t = tncf.time2t(time[:], ttype='bounds')
        np_all_close(t, np.ma.masked_invalid([0, 1, 1, 6, 23]))
        self.assertEqual(t.mask[-1], True)

    def testTime2idx(self):
        from datetime import datetime
        time = np.array([
            datetime.strptime('1970-01-{}+0000'.format(t), '%Y-%m-%d %H:%M%z')
            for t in [
                '01 00:20', '01 01:20', '01 01:40', '01 06:00', '02 01:00'
            ]
        ])
        tncf = self.testncf
        t = tncf.time2idx(time[:])
        np_all_close(t, [0, 1, 2, 6, 23])
        t = tncf.time2idx(
            time, method='bounds', clean='none',
            left=-1, right=999,
        )
        np_all_close(t, [0, 1, 1, 6, 999])
        t = tncf.time2idx(time[:], method='bounds')
        self.assertEqual(True, np.allclose(t, [0, 1, 1, 6, 23]))
        t = tncf.time2idx(time[:], method='bounds', clean='mask', right=np.nan)
        np_all_close(t, np.ma.masked_invalid([0, 1, 1, 6, np.nan]))
        self.assertEqual(t.mask[-1], True)

    @requires_pyproj
    def testXy2ll(self):
        tncf = self.testncf
        crs = tncf.variables['lambert_conformal_conic']
        y0 = crs.false_northing
        x0 = crs.false_easting
        lonmid, latmid = tncf.xy2ll(x0, y0)
        self.assertEqual(True, np.allclose(
            lonmid, crs.longitude_of_central_meridian))
        self.assertEqual(True, np.allclose(
            latmid, crs.latitude_of_projection_origin))
        lon0, lat0 = tncf.xy2ll(0, 0)
        self.assertEqual(True, np.allclose(
            lon0, tncf.variables['longitude_bounds'][0, 0, 0]))
        self.assertEqual(True, np.allclose(
            lat0, tncf.variables['latitude_bounds'][0, 0, 0]))

    @requires_pyproj
    def testIj2ll(self):
        tncf = self.testncf
        lon0, lat0 = tncf.ij2ll(0, 0)
        self.assertEqual(True, np.allclose(
            lon0, tncf.variables['longitude'][0, 0]))
        self.assertEqual(True, np.allclose(
            lat0, tncf.variables['latitude'][0, 0]))

    def testEval(self):
        tncf = self.testncf.copy()
        tncf.eval('O3_PPB = O3 * 1000.', inplace=True, copyall=False)
        o3ppmv = tncf.variables['O3']
        o3ppbv = tncf.variables['O3_PPB']
        self.assertEqual(True, ((o3ppmv * 1000.) == o3ppbv).all())

    def testSetncatts(self):
        tncf = self.testncf.copy()
        tncf.setncatts({'test_new1': 1, 'test_new2': 'five'})
        self.assertEqual(tncf.test_new1, 1)
        self.assertEqual(tncf.test_new2, 'five')

    def testGetncatts(self):
        tncf = self.testncf
        t = tncf.getncatts()
        self.assertEqual(t['str_one'], '1')
        self.assertEqual(t['int_two'], 2)
        self.assertEqual(t['float_threeptfive'], 3.5)

    def testCopy(self):
        tncf = self.testncf
        nncf = tncf.copy(props=True, dimensions=True,
                         variables=True, data=True)
        for pk in tncf.ncattrs():
            self.assertEqual(getattr(tncf, pk), getattr(nncf, pk, None))

        for dk, dv in tncf.dimensions.items():
            dvl = len(dv)
            ndvl = len(nncf.dimensions[dk])
            self.assertEqual(dvl, ndvl)

        for vk, vv in tncf.variables.items():
            nvv = nncf.variables[vk]
            self.assertEqual(True, np.allclose(vv[...], nvv[...]))
            self.assertEqual(vv.dimensions, nvv.dimensions)
            self.assertEqual(vv.dtype.char, nvv.dtype.char)
            for pk in vv.ncattrs():
                pv = getattr(vv, pk)
                npv = getattr(nvv, pk, None)
                testv = pv == npv
                self.assertEqual(True, np.any(testv))

    def testGetTimes(self):
        tncf = self.testncf
        t = tncf.getTimes()
        self.assertEqual(True, (t == self.mytimes).all())

    def testStack(self):
        tncf = self.testncf
        sncf = tncf.stack(tncf, 'TSTEP')
        to3 = tncf.variables['O3'][:]
        no3 = sncf.variables['O3'][:]
        origlen = len(tncf.dimensions['TSTEP'])
        self.assertEqual(origlen * 2, len(sncf.dimensions['TSTEP']))
        self.assertEqual(True, np.allclose(to3, no3[:origlen]))
        self.assertEqual(True, np.allclose(to3, no3[origlen:]))

    def testRemoveSingleton(self):
        tncf = self.testncf
        nncf = tncf.copy()
        nncf.createDimension('test', 1)
        nncf = nncf.removeSingleton(dimkey='test')
        self.assertEqual(set(tncf.dimensions), set(nncf.dimensions))
        nncf.createDimension('test', 1)
        nncf = nncf.removeSingleton(dimkey=None)
        self.assertEqual(set(tncf.dimensions), set(nncf.dimensions))

    def testCreateDimension(self):
        tncf = self.testncf.copy()
        ndims = len(tncf.dimensions)
        olddims = set(tncf.dimensions)
        newdims = olddims.union(['newd'])
        tncf.createDimension('newd', 5)
        self.assertEqual(len(tncf.dimensions), ndims + 1)
        self.assertEqual(set(tncf.dimensions), newdims)

    def testCreateVariable(self):
        tncf = self.testncf.copy()
        nvars = len(tncf.variables)
        var = tncf.createVariable('test', 'f', ('TSTEP',))
        self.assertEqual(len(tncf.variables), nvars + 1)
        self.assertEqual(var.dimensions, ('TSTEP',))
        self.assertEqual(var.dtype.char, 'f')

    def testSave(self):
        import tempfile
        import os
        tncf = self.testncf.copy()
        to3 = tncf.variables['O3'][:]
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmppath = os.path.join(tmpdirname, 'test.nc')
            tncf.save(tmppath).close()
            cncf = pncopen(tmppath)
            np_all_close(to3, cncf.variables['O3'][:])
            sncf = cncf.subset(['O3'])
            np_all_close(to3, sncf.variables['O3'][:])
            os.remove(tmppath)

    def testNcattrs(self):
        tncf = self.testncf
        tncf.ncattrs()

    def testSetncattr(self):
        tncf = self.testncf
        tncf.setncattr('test', 1)
        self.assertEqual(tncf.test, 1)

    def testDelncattr(self):
        tncf = self.testncf.copy()
        tncf.setncattr('test', 1)
        tncf.delncattr('test')
        self.assertEqual(False, hasattr(tncf, 'test'))

    def testInsertDimension(self):
        tncf = self.testncf.copy()
        nncf = tncf.insertDimension(TEST=2, before='LAY')
        checkv = (nncf.variables['O3'][:, 0] == nncf.variables['O3'][:, 1])
        self.assertEqual(True, checkv.all())
        nncf = tncf.insertDimension(TEST=2)
        checkv = (nncf.variables['O3'][0, :] == nncf.variables['O3'][1, :])
        self.assertEqual(True, checkv.all())

    def testInterpDimension(self):
        f1 = PseudoNetCDFFile()
        f1.createDimension('time', 2)
        f1.createDimension('layer', 3)
        f1.createDimension('latitude', 4)
        f1.createDimension('longitude', 5)
        lay = f1.createVariable('layer', 'f', ('layer',))
        lay[:] = np.arange(0, 3)
        simple = f1.createVariable(
            'simple', 'f', ('time', 'layer', 'latitude', 'longitude'))
        simple[0] = np.arange(3 * 4 * 5).reshape(3, 4, 5)
        simple[1] = np.arange(3 * 4 * 5).reshape(3, 4, 5)

        f2 = f1.applyAlongDimensions(layer=lambda x: (x[1:] + x[:-1]) * .5)
        f3 = f1.interpDimension('layer', f2.variables['layer'])
        self.assertEqual(True, np.allclose(
            f2.variables['simple'][:], f3.variables['simple'][:]))
        f4 = PseudoNetCDFFile()
        f4.createDimension('time', 2)
        f4.createDimension('layer', 3)
        f4.createDimension('latitude', 4)
        f4.createDimension('longitude', 5)
        lay = f4.createVariable(
            'layer', 'f', ('time', 'layer', 'latitude', 'longitude'))
        lay[:] = np.arange(0, 3)[None, :, None, None]
        simple = f4.createVariable(
            'simple', 'f', ('time', 'layer', 'latitude', 'longitude'))
        simple[0] = np.arange(3 * 4 * 5).reshape(3, 4, 5)
        simple[1] = np.arange(3 * 4 * 5).reshape(3, 4, 5)

        f5 = f4.applyAlongDimensions(layer=lambda x: (x[1:] + x[:-1]) * .5)
        lay[1] += .25
        f6 = f4.interpDimension('layer', f5.variables['layer'])
        self.assertEqual(True, np.allclose(
            f5.variables['simple'][0], f6.variables['simple'][0]))
        self.assertEqual(False, np.allclose(
            f5.variables['simple'][1], f6.variables['simple'][1]))

    def testNetCDFFileNew(self):
        t = PseudoNetCDFFile.__new__(PseudoNetCDFFile)
        self.assertEqual(t.variables, {})
        self.assertEqual(t.dimensions, {})
        self.assertEqual(t.ncattrs(), ())

    def testNetCDFFileInit(self):
        from numpy import arange
        self._makencf()
        tncf = self.testncf
        self.assertEqual(len(tncf.dimensions['TSTEP']), 24)
        self.assertEqual(len(tncf.dimensions['LAY']), 4)
        self.assertEqual(len(tncf.dimensions['ROW']), 5)
        self.assertEqual(len(tncf.dimensions['COL']), 6)
        self.assertEqual(len(tncf.dimensions['nv']), 4)

        tncf.fish = 2
        setattr(tncf, 'FROG-DOG', 'HAPPY')

        self.assertEqual(set(tncf.variables.keys()), self.myvars)
        o3 = tncf.variables['O3']
        self.assertEqual(
            True, (o3 == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())

        self.assertEqual(o3.typecode(), 'f')

        filedims = list(tncf.dimensions)
        filedims.sort()
        vardims = list(o3.dimensions)
        vardims.sort()
        filedims.remove('nv')
        filedims.remove('tnv')

        self.assertEqual(filedims, vardims)
        from PseudoNetCDF.pncgen import Pseudo2NetCDF
        n = Pseudo2NetCDF().convert(tncf)
        self.assertEqual(set(n.variables.keys()), self.myvars)
        self.assertEqual(
            dict([(k, len(v)) for k, v in n.dimensions.items()]),
            dict([(k, len(v)) for k, v in self.testncf.dimensions.items()])
        )
        self.assertEqual(
            True, (n.variables['O3'][...] == tncf.variables['O3'][...]).all())
        self.assertEqual(n.variables['O3'].units, 'ppbv')
        self.assertEqual(n.fish, 2)
        self.assertEqual(getattr(n, 'FROG-DOG'), 'HAPPY')

    def testNetCDFVariables(self):
        from numpy import arange
        tncf = PseudoNetCDFFile()
        tncf.createDimension('TSTEP', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        def const(*args, **kwds):
            vals = arange(24 * 4 * 5 * 6).reshape((24, 4, 5, 6))
            dims = ('TSTEP', 'LAY', 'ROW', 'COL')
            return PseudoNetCDFVariable(tncf, args[0], 'f', dims, values=vals)
        tncf.variables = PseudoNetCDFVariables(const, ['NO', 'O3'])
        self.assertEqual(True, (tncf.variables['O3'] == arange(
            24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())

    def runTest(self):
        pass
