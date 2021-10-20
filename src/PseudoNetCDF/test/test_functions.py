import unittest


class TestFunc(unittest.TestCase):
    def runTest(self):
        pass

    def setUp(self):
        import PseudoNetCDF as pnc
        import numpy as np
        f = pnc.PseudoNetCDFFile()
        f.createDimension('t', 2)
        f.createDimension('z', 3)
        f.createDimension('y', 4)
        f.createDimension('x', 5)
        var = f.createVariable('test', 'f', ('t', 'z', 'y', 'x'))
        var[:] = np.arange(2 * 3 * 4 * 5).reshape(2, 3, 4, 5)
        self._simplefile = f

    def testPncrename(self):
        from PseudoNetCDF.sci_var import pncrename
        f = self._simplefile.copy()
        nf = pncrename(f, ('d,t,time'))
        chkval = len(f.dimensions['t'])
        outval = len(nf.dimensions['time'])
        self.assertTrue(chkval == outval)
        nf = pncrename(f, ('v,test,hooray'))
        chkval = ['hooray']
        outval = list(nf.variables)
        self.assertTrue(chkval == outval)

    def testManglenames(self):
        from PseudoNetCDF.sci_var import manglenames
        f = self._simplefile.subset([]).copy()
        f.createVariable('arg1$-+', 'f', ('t'))
        nf = manglenames(f)
        outval = list(nf.variables)
        chkval = ['arg1S__add_']
        self.assertTrue(outval == chkval)

    def testRemovesingleton(self):
        from PseudoNetCDF.sci_var import removesingleton
        f = self._simplefile.copy()
        chkval = list(f.dimensions)
        nf = f.copy()
        nf.createDimension('o', 1)
        tmpval = list(nf.dimensions)
        nf2 = removesingleton(nf, 'o')
        outval = list(nf2.dimensions)
        self.assertFalse(tmpval == chkval)
        self.assertTrue(outval == chkval)

    def testSlicedim(self):
        from PseudoNetCDF.sci_var import slice_dim
        import numpy as np
        f = self._simplefile.copy()
        nf = slice_dim(
            slice_dim(
                slice_dim(
                    slice_dim(f, 'x,0,1'),
                    'y,0,1'
                ),
                'z,0,1'
            ),
            't,0,1'
        )
        outval = nf.variables['test'][:]
        chkval = np.array([[[[0]]]])
        self.assertTrue((outval == chkval).all())

    def testReducedim(self):
        from PseudoNetCDF.sci_var import reduce_dim
        import numpy as np
        f = self._simplefile.copy()
        chkval = np.array(f.variables['test'][:]).mean(0, keepdims=True)
        nf = reduce_dim(f, 't,mean')
        outval = nf.variables['test'][:]
        self.assertTrue((outval == chkval).all())

    def testPncbo(self):
        from PseudoNetCDF.sci_var import pncbo
        f = self._simplefile.copy()
        nf = pncbo('+', f, f)
        chkval = f.variables['test'] * 2
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testSeqpncbo(self):
        from PseudoNetCDF.sci_var import seqpncbo
        f = self._simplefile.copy()
        nf = seqpncbo(['+', '+', '+'], [f, f, f, f])
        chkval = f.variables['test'] * 4
        outval = nf[0].variables['test']
        self.assertTrue((outval == chkval).all())

    def testAddattr(self):
        from PseudoNetCDF.sci_var import add_attr
        f = self._simplefile.copy()
        add_attr(f, 'test,global,o,i,1')
        add_attr(f, 'test,test,o,i,1')
        chkval = 1
        outval = f.test
        self.assertTrue(outval == chkval)
        outval = f.variables['test'].test
        self.assertTrue(outval == chkval)
        chkval = 'a'
        add_attr(f, 'test,global,o,c,a')
        add_attr(f, 'test,test,o,c,a')
        outval = f.test
        self.assertTrue(outval == chkval)
        outval = f.variables['test'].test
        self.assertTrue(outval == chkval)
        add_attr(f, 'test,global,d,c,a')
        self.assertFalse(hasattr(f, 'test'))

    def testSplitdim(self):
        from PseudoNetCDF.sci_var import splitdim
        f = self._simplefile.copy()
        tv = f.variables['test']
        chkval = tv.shape[:2] + (2, 2) + tv.shape[-1:]
        nf = splitdim(f, 'y', ('day', 'hour'), (2, 2))
        outval = nf.variables['test'].shape
        self.assertTrue(outval == chkval)

    def testStackfiles(self):
        from PseudoNetCDF.sci_var import stack_files
        import numpy as np
        f = self._simplefile.copy()
        tv = f.variables['test']
        nf = stack_files([f, f, f], 't')
        chkval = np.concatenate([tv] * 3, axis=0)
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testMerge(self):
        from PseudoNetCDF.sci_var import merge
        f = self._simplefile.copy()
        f1 = f.renameVariable('test', 'test1')
        f2 = f.renameVariable('test', 'test2')
        f3 = f.renameVariable('test', 'test3')
        nf = merge([f1, f2, f3])
        chkval = ['test1', 'test2', 'test3']
        outval = sorted(nf.variables)
        self.assertTrue(outval == chkval)

    def testConvolvedim(self):
        from PseudoNetCDF.sci_var import convolve_dim
        f = self._simplefile.copy()
        nf = convolve_dim(f, 't,valid,0.5,0.5')
        chkval = f.variables['test'].mean(0, keepdims=True)
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testPncexpr(self):
        from PseudoNetCDF.sci_var import pncexpr
        f = self._simplefile.copy()
        chkval = f.variables['test'][:] * 2 + 3
        nf = pncexpr('test2 = test * 2 + 3', f)
        outval = nf.variables['test2'][:]
        self.assertTrue((outval == chkval).all())

    def testMaksvals(self):
        from PseudoNetCDF.sci_var import mask_vals
        import numpy as np
        f = self._simplefile.copy()
        chkval = np.ma.masked_greater(f.variables['test'], 5)
        nf = mask_vals(f, 'greater,5')
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testInterpvars(self):
        from PseudoNetCDF.sci_var import interpvars
        import numpy as np
        w = np.array([
            [0.5, 0.5, 0, 0],
            [0, 0.5, 0.5, 0],
            [0, 0, 0.5, 0.5],
        ])
        f = self._simplefile.copy()
        chkval = (
            f.variables['test'][:, :, None, :, :]
            * w[None, None, :, :, None]
        ).sum(3)
        nf = interpvars(f, w, 'y')
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testGetvarpnc(self):
        from PseudoNetCDF.sci_var import getvarpnc
        f = self._simplefile.copy()
        chkval = sorted(f.variables)
        f.copyVariable(f.variables['test'], key='test2')
        tmpval = sorted(f.variables)
        nf = getvarpnc(f, ['test'])
        outval = sorted(nf.variables)
        self.assertFalse(tmpval == chkval)
        self.assertTrue(outval == chkval)

    def testMeshdim(self):
        from PseudoNetCDF.sci_var import mesh_dim
        f = self._simplefile.copy()
        chkval = f.variables['test'][:].reshape(2, 3, 2, 2, 5).mean(3)
        nf = mesh_dim(f, 'y,2,mean')
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testPncfunc(self):
        from PseudoNetCDF.core._functions import pncfunc
        f = self._simplefile.copy()
        chkval = f.variables['test'] * 2 + 7
        nf = pncfunc(lambda x: x * 2 + 7, f)
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())

    def testPncbfunc(self):
        from PseudoNetCDF.core._functions import pncbfunc
        f = self._simplefile.copy()
        chkval = f.variables['test'] * 3
        nf = pncbfunc(lambda x, y: x * 2 + y, f, f)
        outval = nf.variables['test']
        self.assertTrue((outval == chkval).all())


# warnings.warn(str(outval))
"""
def extract_from_file(f, lonlatfs, unique=False, gridded=None, method='nn',
def extract_lonlat(f, lonlat, unique=False, gridded=None, method='nn',
def pncfunc(func, ifile1, coordkeys=None, verbose=0):
def pncbfunc(func, ifile1, ifile2, coordkeys=None, verbose=0):
"""
