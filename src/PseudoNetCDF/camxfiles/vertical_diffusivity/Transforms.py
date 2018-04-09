__all__ = ['vertical_diffusivity_center_time']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx vertical/diffusivity variable transformations
==================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` variable transformations for CAMx
              vertical diffusivity files.  See PseudoNetCDFFile
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from PseudoNetCDF.MetaNetCDF import time_avg_new_unit
from .Memmap import vertical_diffusivity as reg_vertical_diffusivity


class vertical_diffusivity_center_time(time_avg_new_unit):
    __reader__ = reg_vertical_diffusivity
