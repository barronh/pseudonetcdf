from __future__ import print_function
__doc__ = r"""
PseudoNetCDF provides basic interfaces for emulating
NetCDF and manipulating or extending existing NetCDF
like files
"""

__all__ = ['sci_var', 
           'MetaNetCDF',
           'units',
           'ArrayTransforms',
           'camxfiles',
           'test',
           'pncdump',
           'cmaqfiles',
           'racmfiles',
           'geoschemfiles',
           'icarttfiles',
           'toms', 'anyfile', 'conventions']
import sys
import os
import warnings
warn = warnings.warn
def clean_showwarning(message, category, filename, lineno, file=None):
    print('**PNC:%s:%s:%s:\n  %s' % ((filename), lineno, category.__name__, message), file = sys.stderr)
    return
warnings.showwarning = clean_showwarning

from PseudoNetCDF import sci_var
from .sci_var import *
__all__ += sci_var.__all__

from . import camxfiles
from . import cmaqfiles
from . import racmfiles
from . import geoschemfiles
from . import noaafiles
from . import MetaNetCDF
from . import ArrayTransforms
from . import units
from . import icarttfiles
from . import aermodfiles
from . import textfiles
from . import coordutil
from . import test
from ._getreader import anyfile
