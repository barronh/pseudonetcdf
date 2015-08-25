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

__all__ = ['PseudoNetCDFFile', 'PseudoNetCDFDimension', 'PseudoNetCDFVariableConvertUnit', 'PseudoNetCDFFileMemmap', 'PseudoNetCDFVariable', 'PseudoNetCDFMaskedVariable', 'PseudoIOAPIVariable', 'PseudoNetCDFVariables', 'Pseudo2NetCDF', 'reduce_dim', 'slice_dim', 'getvarpnc', 'interpvars', 'extract', 'pncbo', 'seqpncbo', 'pncexpr']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from core._files import PseudoNetCDFFile, PseudoNetCDFFileMemmap, PseudoNetCDFVariables, OrderedDict
from core._dimensions import PseudoNetCDFDimension
from core._variables import PseudoNetCDFVariable, PseudoNetCDFMaskedVariable, PseudoIOAPIVariable
from core._functions import interpvars, extract, mask_vals, slice_dim, reduce_dim, mesh_dim, pncbo, pncexpr, seqpncbo, getvarpnc, add_attr, stack_files, convolve_dim, manglenames, removesingleton, merge
from core._util import get_ncf_object, get_dimension_length
from core._transforms import PseudoNetCDFVariableConvertUnit
from PseudoNetCDF.pncgen import Pseudo2NetCDF
