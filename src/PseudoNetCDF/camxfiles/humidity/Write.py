__all__ = ['ncf2hum']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx humidity  writer
============================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx
              humidity files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from PseudoNetCDF.camxfiles.one3d.Write import ncf2one3d as ncf2humidity

