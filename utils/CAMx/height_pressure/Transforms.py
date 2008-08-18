__all__=['height_pressure_plus', 'height_pressure_center_time_plus', 'height_pressure_center_time']
from numpy import array

from pyPA.utils.MetaNetCDF import add_derived, time_avg_new_unit
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from pyPA.utils.CAMxFiles import height_pressure as reg_height_pressure
from pyPA.utils.ArrayTransforms import CAMxHeightToDepth

class height_pressure_plus(add_derived):
    __childclass__=reg_height_pressure
    __addvars__=['DEPTH']
    def __DEPTH__(self):
        val=CAMxHeightToDepth(self.variables['HGHT'])
        var=PseudoNetCDFVariable(self,'DEPTH','f',('TSTEP','LAY','ROW','COL'),val)
        var.units='m'
        var.long_name='RATE'.ljust(16)
        var.var_desc='RATE'.ljust(16)
        return var

class height_pressure_center_time_plus(time_avg_new_unit):
    __reader__=height_pressure_plus

class height_pressure_center_time(time_avg_new_unit):
    __reader__=reg_height_pressure
    
