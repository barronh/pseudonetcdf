"""
.. xarray Engine

xarray Engine
-------------

PseudoNetCDF can also be used to enable xarray to open non-netcdf files.
This example shows how to open an ICARTT (NASA AMES ffi1001) file, with
xarray. The example file is from the DC-8 during the ARCTAS campaign. Data
includes OH and HO2 volume mixing ratios (VMR) from 2008-06-28T14:22:18Z to
2008-06-28T21:22:36Z with 7s samples.

.. code-block:: python
"""
if True:
  import xarray as xr
  import urllib.request
  import os

  hoxpath = 'HOX_DC8_2008-06-26_r2008-12-21.ict'
  # Get the file if -- url may change.
  if not os.path.exists(hoxpath):
    url = 'http://www-air.larc.nasa.gov/cgi-bin/enzFile?c16141B08DF7F1ACFBAD5'
    url += 'C83F9313E20C792f7075622d6169722f4152435441532f4443385f41495243524'
    url += '146542f4252554e452e57494c4c49414d2f484f785f4443385f32303038303632'
    url += '365f52312e696374'
    urllib.request.urlretrieve(url, hoxpath)

  opts = dict(engine='pseudonetcdf', backend_kwargs=dict(format='ffi1001'))
  f = xr.open_dataset(hoxpath, **opts)
  print(f)
  # <xarray.Dataset>
  # Dimensions:    (POINTS: 1188)
  # Dimensions without coordinates: POINTS
  # Data variables:
  #     Start_UTC  (POINTS) timedelta64[ns] ...
  #     Stop_UTC   (POINTS) timedelta64[ns] ...
  #     Mid_UTC    (POINTS) timedelta64[ns] ...
  #     OH_pptv    (POINTS) float64 ...
  #     HO2_pptv   (POINTS) float64 ...
  # Attributes: (12/31)
  #     fmt:                              1001
  #     n_header_lines:                   36
  #     PI_NAME:                          Brune, William
  #     ORGANIZATION_NAME:                Penn State University
  #     SOURCE_DESCRIPTION:               ATHOS - OH and HO2 concentrations using...
  #     MISSION_NAME:                     NASA ARCTAS Mission 2008
  #     ...                               ...
  #     PROJECT_INFO:                     NASA ARCTAS Mission 1 April-15 July 200...
  #     STIPULATIONS_ON_USE:              Use of these data requires prior approv...
  #     OTHER_COMMENTS:                   These Data are PRELIMINARY; contact the...
  #     REVISION:                         R1
  #     R1:                               Final Data - StartUTC and StopUTC of ea...
  #     TFLAG:                            Start_UTC