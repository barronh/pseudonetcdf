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
           'toms',
           'anyfile',
           'pncopen',
           'pncwrite',
           'conventions',
           'woudcfiles']
import sys
import os
from PseudoNetCDF.pncwarn import warn

def makequite():
    global _quiet
    _quiet = True
from PseudoNetCDF import sci_var
from .sci_var import *
__all__ += sci_var.__all__

from . import camxfiles
from . import cmaqfiles
from . import racmfiles
from . import geoschemfiles
from . import noaafiles
from . import woudcfiles
from . import epafiles
from . import MetaNetCDF
from . import ArrayTransforms
from . import units
from . import icarttfiles
from . import aermodfiles
from . import textfiles
from . import ceilometerfiles
from . import coordutil
from . import test
from ._getreader import anyfile, pncopen
from ._getwriter import pncwrite
from .pncparse import PNC, pnc
