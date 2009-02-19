__all__=['height_pressure_plus', 'height_pressure_center_time_plus', 'height_pressure_center_time']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx height/pressure variable transformations
=============================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` variable transformations
              for CAMx height/pressure files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""


from numpy import array

from PseudoNetCDF.MetaNetCDF import add_derived, time_avg_new_unit
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from PseudoNetCDF.camxfiles.height_pressure.Memmap import height_pressure as reg_height_pressure
from PseudoNetCDF.ArrayTransforms import CAMxHeightToDepth

class height_pressure_plus(add_derived):
    __childclass__=reg_height_pressure
    __addvars__=['DEPTH']
    def __DEPTH__(self):
        val=CAMxHeightToDepth(self.variables['HGHT'])
        var=PseudoNetCDFVariable(self,'DEPTH','f',('TSTEP','LAY','ROW','COL'),values=val)
        var.units='m'
        var.long_name='RATE'.ljust(16)
        var.var_desc='RATE'.ljust(16)
        return var

class height_pressure_center_time_plus(time_avg_new_unit):
    __reader__=height_pressure_plus

class height_pressure_center_time(time_avg_new_unit):
    __reader__=reg_height_pressure
    
