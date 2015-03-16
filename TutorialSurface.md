# sfc\_pres\_temp #

This tutorial mimics the sfc\_pres\_temp tutorial on the NetCDF page.

```
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from PseudoNetCDF.pncdump import pncdump
from numpy import arange

sfc_pres_temp_file = PseudoNetCDFFile()

sfc_pres_temp_file.createDimension('latitude', 6)
sfc_pres_temp_file.createDimension('longitude', 12)

# Create the latitude variable
lat = sfc_pres_temp_file.createVariable('latitude', 'f', ('latitude',))

# Add the units property
lat.units = 'degrees_north'

# Add data
lat[:] = [25, 30, 35, 40, 45,50]

# Create the longitude variable
lon = sfc_pres_temp_file.createVariable('longitude', 'f', ('longitude',))

# Add the units property
lon.units = 'degrees_east'

# Add the data
lon[:] = [-125, -120, -115, -110, -105, -100, -95, -90, -85, -80, -75, -70]

# Create the pressure variable with dimensions (latitude, longitude)
press = sfc_pres_temp_file.createVariable('pressure', 'f', ('latitude', 'longitude'))

# Add the units property
press.units = 'hPa'

#Create the pressure values
press_vals = arange(900,972).reshape(12,6).swapaxes(0,1)

# Assign the pressure values to the variable
press[:] = press_vals


# Create the temperature variable with dimensions (latitude, longitude)
temp = sfc_pres_temp_file.createVariable('temperature', 'f', ('latitude', 'longitude'))

# Add the units property
temp.units = 'celsius'

#Create the temperature values
temp_vals = arange(9,27.,.25).reshape(12,6).swapaxes(0,1)

# Assign the temperature values to the variable
temp[:] = temp_vals

pncdump(sfc_pres_temp_file, name = 'sfc_pres_temp')
```