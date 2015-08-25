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
           'toms', 'getfile', 'conventions']
import sys
import os
import warnings
warn = warnings.warn
def clean_showwarning(message, category, filename, lineno, file=None):
    print >> sys.stderr, '**PNC:%s:%s:%s:\n  %s' % ((filename), lineno, category.__name__, message)
    return
warnings.showwarning = clean_showwarning

import sci_var
from sci_var import *
__all__ += sci_var.__all__

import camxfiles
import cmaqfiles
import racmfiles
import geoschemfiles
import noaafiles
import MetaNetCDF
import ArrayTransforms
import units
import icarttfiles
import aermodfiles
import textfiles
import coordutil
import test
from _getreader import anyfile
