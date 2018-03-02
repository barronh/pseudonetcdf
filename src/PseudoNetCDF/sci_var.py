__doc__ = r"""
Scientific data has dimensions that have physical meaning and values
only have meaning in the context of their units.  This module implements
numpy arrays that are aware of their dimensions trying to vaguely adhere
to the Common Data Model from Unitdata at UCAR.

Each variable has as a property its dimensions names (dimensions).  Further,
each dimension name exists as a property and contains a one dimensional array
of values associated with that dimension.

For the purposes of ease of use, the standard properties of netCDF files
are attached and the arrays implement the Scientific.IO.NetCDF.NetCDFVariable
interfaces.
"""

__all__ = ['PseudoNetCDFFile', 'PseudoNetCDFDimension',
           'PseudoNetCDFVariableConvertUnit',
           'PseudoNetCDFVariable',
           'PseudoNetCDFMaskedVariable',
           'PseudoIOAPIVariable',
           'PseudoNetCDFVariables',
           'Pseudo2NetCDF',
           'reduce_dim',
           'slice_dim',
           'getvarpnc',
           'interpvars',
           'extract',
           'pncbo',
           'seqpncbo',
           'pncexpr',
           'WrapPNC',
           'manglenames',
           'merge',
           'pncrename',
           'removesingleton',
           'splitdim',
           'mask_vals',
           'add_attr',
           'stack_files',
           'extract_from_file',
           'mesh_dim',
           'convolve_dim',
           'get_ncf_object',
           'get_dimension_length']


from .core._files import PseudoNetCDFFile, PseudoNetCDFVariables
from .core._wrapnc import WrapPNC
from .core._dimensions import PseudoNetCDFDimension
from .core._variables import PseudoNetCDFVariable, PseudoNetCDFMaskedVariable
from .core._variables import PseudoIOAPIVariable
from .core._functions import interpvars, extract, mask_vals, slice_dim
from .core._functions import reduce_dim, mesh_dim, pncbo, pncexpr, seqpncbo
from .core._functions import getvarpnc, add_attr, stack_files, convolve_dim
from .core._functions import manglenames, removesingleton, merge
from .core._functions import extract_from_file, pncrename, splitdim
from .core._util import get_ncf_object, get_dimension_length
from .core._transforms import PseudoNetCDFVariableConvertUnit
from PseudoNetCDF.pncgen import Pseudo2NetCDF
