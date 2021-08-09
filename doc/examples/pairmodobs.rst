.. Pair Model with AQS

Pair Model Data with AQS
~~~~~~~~~~~~~~~~~~~~~~~~

It is common to evaluate Air Quality Models with AQS, so the first example
uses the TestCase from CMAQ to do just that. First, download the `testcase
for CMAQ <https://www.epa.gov/cmaq/cmaq-inputs-and-test-case-data>`_ and
the hourly data from AQS `Pre-generated Hourly Table 
<https://aqs.epa.gov/aqsweb/airdata/download_files.html>`_ for 
ozone. Unzip both to a working directory.

.. code-block:: python

  import PseudoNetCDF as pnc
  import numpy as np
  import pandas as pd

  aqpath = 'CCTM_ACONC_v52_cb6r3_intel17.0_SE52BENCH_20110701.nc'
  aqf = pnc.pncopen(aqpath, format='ioapi').subsetVariables(['O3'])

  aqspath = 'hourly_44201_2011.csv'
  aqsdata = pd.read_csv(aqspath, parse_dates = [["Date GMT", "Time GMT"]])

  lat = aqsdata.Latitude.values
  lon = aqsdata.Longitude.values
  
  # Convert time to hourly steps for each observation
  starttime = aqf.getTimes()[0].replace(tzinfo=None)
  time = aqsdata['Date GMT_Time GMT'].dt.to_pydatetime()
  dt = np.array([(t - starttime).total_seconds() // 3600 for t in time], 'i')
  aqsdata['T'] = dt
  
  # Convert lon/lat to indices and add a surface index
  aqsdata['I'], aqsdata['J'] = aqf.ll2ij(lon, lat)
  aqsdata['K'] = aqsdata['I'] * 0
  
  # Only include data in time/space domain
  nx = aqf.NCOLS
  ny = aqf.NROWS
  nt = len(aqf.dimensions['TSTEP'])
  aqssubset = aqsdata.query(
      f'I > 0 and I < {nx} and J > 0 and J < {ny} and T > 0 and T < {nt}'
  ).copy()

  # Index model ozone and add to obs DataFrame
  i = aqssubset['I'].values
  j = aqssubset['J'].values
  k = aqssubset['K'].values
  t = aqssubset['T'].values
  aqssubset['ModOzone'] = aqf.variables['O3'][:][t, k, j, i]

  # Print Mean and Correlation
  print(aqssubset.filter(['Sample Measurement', 'ModOzone']).mean())
  print(aqssubset.filter(['Sample Measurement', 'ModOzone']).corr())
  
  # Make a scatter plot
  ax = aqssubset.plot.scatter(x='Sample Measurement', y='ModOzone')
  ax.figure.savefig('Ozone.png')

  # Output printed from above
  # Sample Measurement    0.041571
  # ModOzone              0.044702
  # dtype: float64
  #                     Sample Measurement  ModOzone
  # Sample Measurement            1.000000  0.816638
  # ModOzone                      0.816638  1.000000
  # Scatter plot figure not shown

This example will work just as well with a CAMx file and the uamiv reader, or
a GEOS-Chem file and the bpch reader.
