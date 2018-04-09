__all__ = ['ncf2vertical_diffusivity']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx vertical diffusivity  writer
=================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx vertical
              diffusivity files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from PseudoNetCDF.camxfiles.one3d.Write import ncf2one3d
from PseudoNetCDF._getwriter import registerwriter

ncf2vertical_diffusivity = ncf2one3d
registerwriter('camxfiles.vertical_diffusivity', ncf2vertical_diffusivity)
registerwriter('vertical_diffusivity', ncf2vertical_diffusivity)
