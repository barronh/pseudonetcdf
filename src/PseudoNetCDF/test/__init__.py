import importlib
import unittest
from distutils.version import LooseVersion


def _importorskip(modname, minversion=None):
    try:
        mod = importlib.import_module(modname)
        has = True
        if minversion is not None:
            if LooseVersion(mod.__version__) < LooseVersion(minversion):
                raise ImportError('Minimum version not satisfied')
    except ImportError:
        has = False
    # TODO: use pytest.skipif instead of unittest.skipUnless
    # Using `unittest.skipUnless` is a temporary workaround for pytest#568,
    # wherein class decorators stain inherited classes.
    # xref: xarray#1531, implemented in xarray #1557.
    func = unittest.skipUnless(has, reason='requires {}'.format(modname))
    return has, func


has_pyproj, requires_pyproj = _importorskip('pyproj')
has_matplotlib, requires_matplotlib = _importorskip('matplotlib')
has_basemap, requires_basemap = _importorskip('mpl_toolkits.basemap')


def compare_files(f1, f2):
    from numpy.testing import assert_array_equal
    assert (set(f1.variables) == set(f2.variables))
    assert (set(f1.dimensions) == set(f2.dimensions))
    attrs1 = f1.getncatts()
    attrs2 = f2.getncatts()
    assert (set(attrs1) == set(attrs2))
    for key, val1 in f1.dimensions.items():
        val2 = f2.dimensions[key]
        assert_array_equal(len(val1), len(val2))
        assert_array_equal(val1.isunlimited(), val2.isunlimited())

    for key, val1 in attrs1.items():
        val2 = attrs2[key]
        assert_array_equal(val1, val2)

    for key, val1 in f1.variables.items():
        val2 = f2.variables[key]
        assert_array_equal(val1[...], val2[...])
