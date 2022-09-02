from __future__ import print_function
__doc__ = r"""
PseudoNetCDF provides basic interfaces for emulating
NetCDF and manipulating or extending existing NetCDF
like files
"""

__all__ = ['aermodfiles',
           'anyfile',
           'ArrayTransforms',
           'camxfiles',
           'ceilometerfiles',
           'cmaqfiles',
           'conventions',
           'coordutil',
           'epafiles',
           'pncfile',
           'pncvar',
           'geoschemfiles',
           'icarttfiles',
           'MetaNetCDF',
           'noaafiles',
           'pnc',
           'PNC',
           'pncdump',
           'pncopen',
           'pncmfopen',
           'getreader',
           'getreaderdict',
           'pncwrite',
           'PseudoNetCDFFile',
           'PseudoNetCDFVariable',
           'PseudoNetCDFVariables',
           'racmfiles',
           'sci_var',
           'test',
           'textfiles',
           'toms',
           'units',
           'warn',
           'woudcfiles',
           'wrffiles',
           '__version__']

from .pncwarn import warn
from .version import version as __version__

from .core import PseudoNetCDFFile
from .core import PseudoNetCDFVariable, PseudoNetCDFVariables

from . import sci_var
from . import units
from . import test
from . import coordutil
from . import MetaNetCDF
from . import ArrayTransforms

from . import aermodfiles
from . import camxfiles
from . import ceilometerfiles
from . import cmaqfiles
from . import epafiles
from . import geoschemfiles
from . import icarttfiles
from . import noaafiles
from . import racmfiles
from . import textfiles
from . import toms
from . import woudcfiles
from . import wrffiles

from ._getreader import anyfile, pncopen, pncmfopen, getreader, getreaderdict
from ._getwriter import pncwrite
from .pncparse import PNC, pnc

pncfile = PseudoNetCDFFile
pncvar = PseudoNetCDFVariable

for k in sci_var.__all__:
    if k not in __all__:
        __all__.append(k)
    globals()[k] = getattr(sci_var, k)


def makequite():
    global _quiet
    _quiet = True
