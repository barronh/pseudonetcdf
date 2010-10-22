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
           'toms']
import sci_var
from sci_var import *
__all__ += sci_var.__all__

import MetaNetCDF
import ArrayTransforms
import units
import camxfiles
import cmaqfiles
import racmfiles
import icarttfiles
import test
