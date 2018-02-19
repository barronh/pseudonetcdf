import importlib
import unittest

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
has_basemap, requires_basemap = _importorskip('mpl_toolkits.basemap')
