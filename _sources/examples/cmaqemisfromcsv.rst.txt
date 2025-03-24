.. CMAQ Emissions from CSV

Make CMAQ Emissions from a CSV File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make a CMAQ emissions input for NOx from a CSV file with latitude,
longitude, and hour-of-day coordinates.  The code below will make a file from
a `GRIDDESC` file. Then, add a variable for NO and NO2. Finally, populate the
file and save it as an IOAPI-like file.

.. code-block:: python

  import PseudoNetCDF as pnc
  import pandas as pd
  import numpy as np

  # create a dummy input file
  lat = (np.sin(np.linspace(0, 2*np.pi, 100)) * 7 + 40).round(2)
  lon = np.linspace(-120, -80, 100)
  hour = np.floor(np.linspace(0, 23.99, 100)).astype('i')
  NO = np.random.lognormal(0.09, sigma=0.1, size=100)
  NO2 = np.random.lognormal(0.01, sigma=0.01, size=100)
  data = pd.DataFrame.from_dict(dict(
      lat=lat, lon=lon, hour=hour, NO=NO, NO2=NO2
  )).to_csv('emis.csv', index=False)
  
  # Either open GRIDDESC from disk or create one in memory
  gdpath = '/path/to/GRIDDESC'
  gdpath = """' '
  'LamCon_40N_97W'\n 2 33.000 45.000 -97.000 -97.000 40.000
  ' '
  '12US1'
  'LamCon_40N_97W' -2556000.0 -1728000.0 12000.0 12000.0 459 299 1
  ' '
  """
  
  # Create a Template
  gf = pnc.pncopen(gdpath, format='griddesc')
  gf.SDATE = 2016001
  gf.TSTEP = 10000
  gf.updatetflag(overwrite=True)
  tmpf = gf.sliceDimensions(TSTEP=[0]*25)

  #  Read in the emissions
  data = pd.read_csv('emis.csv')

  # Add I,J and time coordinate
  data['I'], data['J'] = tmpf.ll2ij(data.lon.values, data.lat.values)
  cellhourtotal = data.groupby(['I', 'J', 'hour'], as_index=False).sum()
  # Track emission keys
  metakeys = ('lat', 'lon', 'hour')
  emiskeys = [k for k in data.columns if k not in metakeys]

  h = cellhourtotal.hour
  i = cellhourtotal.I
  j = cellhourtotal.J

  for emiskey in emiskeys:
      evar = tmpf.createVariable(emiskey, 'f', ('TSTEP', 'LAY', 'ROW', 'COL'))
      evar.setncatts(dict(
        units='moles/s', long_name=emiskey, var_desc=emiskey
      ))
      evar[h, h*0, j, i] = cellhourtotal[emiskey].values

  # Get rid of initial DUMMY variable
  del tmpf.variables['DUMMY']

  # Update TFLAG to be consistent with variables
  tmpf.updatetflag(tstep=10000, overwrite=True)

  # Remove VAR-LIST so that it can be inferred
  delattr(tmpf, 'VAR-LIST')
  tmpf.updatemeta()

  # Save out
  tmpf.save('emis.nc', format='NETCDF3_CLASSIC')

