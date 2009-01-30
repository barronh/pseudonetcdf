__all__=['avrg_conc', \
         'cloud_rain', \
         'depn', \
         'one3d', \
         'one3d_mem', \
         'gridded_emissions', \
         'height_pressure', \
         'height_pressure_mem', \
         'humidity', \
         'humidity_mem', \
         'inst_conc', \
         'ipr', \
         'ipr2ncf', \
         'ipr_mem', \
         'irr', \
         'irr2ncf', \
         'irr_mem', \
         'landuse', \
         'lo_emiss', \
         'point_source', \
         'point_source_mem', \
         'temperature', \
         'temperature_mem', \
         'uamiv', \
         'uamiv_mem', \
         'vertical_diffusivity', \
         'vertical_diffusivity_mem', \
         'wind', \
         'wind_mem', \
         'write_emissions', \
         'write_emissions_ncf', \
         'write_hgtprss', \
         'write_point', \
         'write_wind']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import unittest

from PseudoNetCDF.camxfiles.cloud_rain.Memmap import cloud_rain
from PseudoNetCDF.camxfiles.height_pressure.Memmap import height_pressure
from PseudoNetCDF.camxfiles.height_pressure.Read import height_pressure as height_pressure_mem
from PseudoNetCDF.camxfiles.height_pressure.Write import write_hgtprss
from PseudoNetCDF.camxfiles.humidity.Memmap import humidity
from PseudoNetCDF.camxfiles.humidity.Read import humidity as humidity_mem
from PseudoNetCDF.camxfiles.ipr.Memmap import ipr
from PseudoNetCDF.camxfiles.ipr.Read import ipr as ipr_mem
from PseudoNetCDF.camxfiles.irr.Memmap import irr
from PseudoNetCDF.camxfiles.irr.Read import irr as irr_mem
from PseudoNetCDF.camxfiles.landuse.Memmap import landuse
from PseudoNetCDF.camxfiles.one3d.Memmap import one3d
from PseudoNetCDF.camxfiles.one3d.Read import one3d as one3d_mem
from PseudoNetCDF.camxfiles.point_source.Memmap import point_source
from PseudoNetCDF.camxfiles.point_source.Read import point_source as point_source_mem
from PseudoNetCDF.camxfiles.point_source.Write import write_point
from PseudoNetCDF.camxfiles.temperature.Memmap import temperature
from PseudoNetCDF.camxfiles.temperature.Read import temperature as temperature_mem
from PseudoNetCDF.camxfiles.uamiv.Memmap import uamiv
from PseudoNetCDF.camxfiles.uamiv.Read import uamiv as uamiv_mem
from PseudoNetCDF.camxfiles.uamiv.Write import write_emissions_ncf,write_emissions
from PseudoNetCDF.camxfiles.vertical_diffusivity.Memmap import vertical_diffusivity
from PseudoNetCDF.camxfiles.vertical_diffusivity.Read import vertical_diffusivity as vertical_diffusivity_mem
from PseudoNetCDF.camxfiles.wind.Memmap import wind
from PseudoNetCDF.camxfiles.wind.Read import wind as wind_mem
from PseudoNetCDF.camxfiles.wind.Write import write_wind

lo_emiss=uamiv
avrg_conc=uamiv
inst_conc=uamiv
depn=uamiv
gridded_emissions=uamiv

from CAMxFileConverter import irr2ncf,ipr2ncf

__doc__="""
All file interfaces should implement a Scientific.IO.NetCDF.NetCDFFile like interface.

Any other methods are purely for internal use and should not be expected to perform
consistently from one revision to the next.  The methods are not class only variables for
backward compatibility that will not be supported.
"""

if __name__ == '__main__':
    unittest.main()
