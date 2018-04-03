__all__ = ['humidity_center_time']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx humidity variable transformations
==================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` variable transformations
              for CAMx humidity files.  See PseudoNetCDFFile
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from PseudoNetCDF.MetaNetCDF import time_avg_new_unit
from .Memmap import humidity as reg_humidity


class humidity_center_time(time_avg_new_unit):
    __reader__ = reg_humidity
