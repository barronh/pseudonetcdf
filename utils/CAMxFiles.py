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
         'unittest', \
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
from CAMx.cloud_rain.Memmap import cloud_rain
from CAMx.height_pressure.Memmap import height_pressure
from CAMx.height_pressure.Read import height_pressure as height_pressure_mem
from CAMx.height_pressure.Write import write_hgtprss
from CAMx.humidity.Memmap import humidity
from CAMx.humidity.Read import humidity as humidity_mem
from CAMx.ipr.Memmap import ipr
from CAMx.ipr.Read import ipr as ipr_mem
from CAMx.irr.Memmap import irr
from CAMx.irr.Read import irr as irr_mem
from CAMx.landuse.Memmap import landuse
from CAMx.one3d.Memmap import one3d
from CAMx.one3d.Read import one3d as one3d_mem
from CAMx.point_source.Memmap import point_source
from CAMx.point_source.Read import point_source as point_source_mem
from CAMx.point_source.Write import write_point
from CAMx.temperature.Memmap import temperature
from CAMx.temperature.Read import temperature as temperature_mem
from CAMx.uamiv.Memmap import uamiv
from CAMx.uamiv.Read import uamiv as uamiv_mem
from CAMx.uamiv.Write import write_emissions_ncf,write_emissions
from CAMx.vertical_diffusivity.Memmap import vertical_diffusivity
from CAMx.vertical_diffusivity.Read import vertical_diffusivity as vertical_diffusivity_mem
from CAMx.wind.Memmap import wind
from CAMx.wind.Read import wind as wind_mem
from CAMx.wind.Write import write_wind

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
