__all__ = ['temperature_center_time']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx temperature variable transformations
=========================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` variable transformations
              for CAMx temperature files.  See pnc.sci_var.PseudoNetCDFFile
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from PseudoNetCDF.MetaNetCDF import time_avg_new_unit
from .Memmap import temperature as reg_temperature


class temperature_center_time(time_avg_new_unit):
    __reader__ = reg_temperature
