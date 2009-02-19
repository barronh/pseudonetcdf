__all__=['temperature_center_time']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx temperature variable transformations
=========================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` variable transformations
              for CAMx temperature files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from numpy import array

from PseudoNetCDF.MetaNetCDF import add_derived, time_avg_new_unit
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from Memmap import temperature as reg_temperature

class temperature_center_time(time_avg_new_unit):
    __reader__=reg_temperature
    
