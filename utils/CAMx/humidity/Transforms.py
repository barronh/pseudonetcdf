__all__=['humidity_center_time']
from numpy import array

from pyPA.utils.MetaNetCDF import add_derived, time_avg_new_unit
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from pyPA.utils.CAMxFiles import humidity as reg_humidity

class humidity_center_time(time_avg_new_unit):
    __reader__=reg_humidity
