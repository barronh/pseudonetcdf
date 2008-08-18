__all__=['vertical_diffusivity_center_time']
from numpy import array

from pyPA.utils.MetaNetCDF import add_derived, time_avg_new_unit
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from pyPA.utils.CAMxFiles import vertical_diffusivity as reg_vertical_diffusivity

class vertical_diffusivity_center_time(time_avg_new_unit):
    __reader__=reg_vertical_diffusivity
    
