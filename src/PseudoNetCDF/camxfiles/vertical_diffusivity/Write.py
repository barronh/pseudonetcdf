__all__ = ['ncf2kv']
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

from PseudoNetCDF.camxfiles.one3d.Write import ncf2one3d as ncf2vertical_diffusivity
