.. Examples

Example Use Cases
-----------------

PseudoNetCDF gets used for hundreds of different activities. So far, there are
two total usage examples used to illustrate its power.


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
  
Process Air Quality Forecast Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One user asked for code to process forecasting data where they wanted to

#. Subset times from two daily CMAQ files
#. Concatenate subsets
#. Combine with similarly subset data from WRF
#. Derive a variable
#. Make plots

.. code-block:: python

  import PseudoNetCDF as pnc

  # Read in files and slice times and layers
  aqpath1 = 'CCTM_ACONC_20141225'
  aqpath2 = 'CCTM_ACONC_20141226'
  aqf1 = pnc.pncopen(aqpath1, format='ioapi')\
                .subsetVariables(['O3'])\
                .sliceDimensions(TSTEP=slice(-5, None), LAY=0)
  aqf2 = pnc.pncopen(aqpath2, format='ioapi')\
                .subsetVariables(['O3'])\
                .sliceDimensions(TSTEP=slice(None,19), LAY=0)
  
  # Stack on TSTEP
  # also useful with lists aqf = aqfs[0].stack(aqfs[1:], 'TSTEP')
  aqf = aqf1.stack(aqf2, 'TSTEP')

  # Same with WRF, but fitting to AQ domain. Not actually the same, but doing it any way
  wrfpath = 'wrfout_d01_2014-12-24_12:00:00'
  wrff = pnc.pncopen(wrfpath)\
                .subsetVariables(['T2'])\
                .sliceDimensions(
                    Time=slice(31, 55),
                    bottom_top=0,
                    west_east=slice(100, 287),
                    south_north=slice(100, 287)
                )


  aqf.copyVariable(
    wrff.variables['T2'],
    key='T2', dimensions=('TSTEP', 'LAY', 'ROW', 'COL')
  )

  aqf.eval('O3overT2 = O3/T2', inplace=True)
  ax = aqf.plot('O3overT2')
  ax.figure.savefig('O3overT2.png')
