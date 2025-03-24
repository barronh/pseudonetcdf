.. ICARTT Plots and Conversion

ICARTT Plots and Conversion
---------------------------

This example shows how to open an ICARTT (NASA AMES ffi1001) file, resample
it in time, and then save it out as a NetCDF file. The example file is from
the DC-8 during the ARCTAS campaign. Data includes OH and HO2 volume mixing
ratios (VMR) from 2008-06-28T14:22:18Z to 2008-06-28T21:22:36Z with 7s samples.

.. code-block:: python

  import PseudoNetCDF as pnc
  import matplotlib.pyplot as plt
  import urllib.request
  import os
  import pandas as pd
  
  hoxpath = 'HOX_DC8_2008-06-26_r2008-12-21.ict'
  # Get the file if -- url may change.
  if not os.path.exists(hoxpath):
    url = 'http://www-air.larc.nasa.gov/cgi-bin/enzFile?c16141B08DF7F1ACFBAD5'
    url += 'C83F9313E20C792f7075622d6169722f4152435441532f4443385f41495243524'
    url += '146542f4252554e452e57494c4c49414d2f484f785f4443385f32303038303632'
    url += '365f52312e696374'
    urllib.request.urlretrieve(url, hoxpath)

  # Open the file and plot directly
  f = pnc.pncopen(hoxpath, format='ffi1001')
  x = f.variables['Start_UTC'][:]
  fig, ax = plt.subplots()
  ax.scatter(x, f.variables['OH_pptv'][:] * 100, label='OH VMR x1e14')
  ax.scatter(x, f.variables['HO2_pptv'][:], label='HO2 VMR x1e12')
  ax.legend()
  fig.savefig('icartt_hox.png')
  
  tkeys = ['Start_UTC', 'Mid_UTC', 'Stop_UTC']
  vkeys = ['OH_pptv', 'HO2_pptv']
  ckeys = ['OH_count', 'HO2_count']
  keys = tkeys + vkeys

  # Create a 300s merge file
  freq = '300s'
  # convert PseudoNetCDF to DataFrame
  df = pd.DataFrame({k: f.variables[k] for k in keys})
  # add time element
  refdate = pd.to_datetime('2008-06-28T00Z')
  df['time'] = refdate + pd.to_timedelta(df['Mid_UTC'], unit='s')
  # group by time at specified frequency
  mdg = df.groupby(pd.Grouper(key='time', freq=freq))
  ops = {k: (k, 'mean') for k in keys}
  ops['OH_count_1'] = ('OH_pptv', 'count')
  ops['HO2_count_1'] = ('HO2_pptv', 'count')
  ops['Mid_count_1'] = ('Mid_UTC', 'count')
  mdf = mdg.agg(**ops)

  # Convert to a PseudoNetCDFFile to a dictionary of arrays for convenience
  mf = pnc.PseudoNetCDFFile()
  mf.createDimension('POINTS', mdf.shape[0])
  for k, s in mdf.items():
    vk, _, unit = k.rpartition('_')
    print(vk, s)
    mf.createVariable(vk, 'd', ('POINTS',), units=unit)[:] = s.values
  # Save the file as a NetCDF file
  mf.save(f'{hoxpath}_{freq}.nc')
  
  # Open the netcdf and plot
  nf = pnc.pncopen(f'{hoxpath}_{freq}.nc', format='netcdf')
  x = nf.variables['Start'][:]
  fig, ax = plt.subplots()
  ax.scatter(x, nf.variables['OH'][:] * 100, label='OH VMR x1e14')
  ax.scatter(x, nf.variables['HO2'][:], label='HO2 VMR x1e12')
  ax.legend()
  fig.savefig(f'{freq}_hox.png')
