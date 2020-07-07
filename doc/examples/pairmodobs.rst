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
  import pandas as pd

  aqpath = 'CCTM_ACONC_v52_cb6r3_intel17.0_SE52BENCH_20110701.nc'
  aqf = pnc.pncopen(aqpath, format='ioapi')
  
  aqspath = 'hourly_44201_2011.csv'
  aqsdata = pd.read_csv(aqspath, parse_dates = [["Date GMT","Time GMT"]])
  
  lat = aqsdata.Latitude.values
  lon = aqsdata.Longitude.values
  time = aqsdata['Date GMT_Time GTM'].values
  i, j = aqf.ll2ij(lon, lat)
  k = i * 0
  t = aqf.time2t(time)
  aqsdata['ModOzone'] = aqsdata.sliceDimensions(TSTEP=t, LAY=k, ROW=j, COL=i)
  
  ax = aqsdata.plot(x='Ozone', y='ModOzone')
  ax.figure.savefig('Ozone.png')
  
This example will work just as well with a CAMx file and the uamiv reader, or
a GEOS-Chem file and the bpch reader.
