.. Quick

Quick Start
-----------

This section provides a quick overview of how to work with PseudoNetCDF. Each
quick example is cumulative. If you don't understand something, look at the 
previous sections.

Import PseudoNetCDF
~~~~~~~~~~~~~~~~~~~

First, PseudoNetCDF should be imported as `pnc` for brevity.

.. code-block:: python

  import PseudoNetCDF as pnc

Open a File
~~~~~~~~~~~

:func:`~PseudoNetCDF.pncopen` helps PseudoNetCDF open many types of files, but
is most well tested with Air Quality related files including models (e.g.,
CAMx, CMAQ, GEOS-Chem) and observations (e.g., woudc sondes, AQS, aircraft
icartt). These examples are from publicly available data.

.. code-block:: python

  # from CAMx Test Case v6
  camxpath1 = 'CAMx.v6.40.midwest.36.12.noMPI.20020603.avrg.grd02'
  camxf = pnc.pncopen(camxpath1, format='uamiv')

.. code-block:: python

  # from CMAQ Test Case v5.2
  cmaqpath1 = 'CCTM_ACONC_v52_cb6r3_intel17.0_SE52BENCH_20110701.nc'
  cmaqf = pnc.pncopen(cmaqpath1, format='ioapi')

.. code-block:: python

  # from the GEOS-Chem benchmark v12
  gcpath1 = 'ctm.bpch'
  gcf = pnc.pncopen(gcpath1, format='bpch')

.. code-block:: python

  # from the INTEX-NA NASA campaign
  ictpath = 'HOX_DC8_20040626_R0.ict'
  ictf = pnc.pncopen(ictpath, format='ffi1001')

.. code-block:: python

  # from the woudc.org
  woudcpath = 'bu20170609.b18.csv'
  wdcf = pnc.pncopen(woudcpath, format='woudcsonde')

Don't know the format of your file... leave out the format keyword and
PseudoNetCDF will try to figure it out for you.

Get Help
~~~~~~~~

Most objects have are documented, so use help if you want to know more. For
example, if you wanted to learn more about opening a file
:func:`pnc.pncopen`, use :func:help

.. code-block:: python

  help(pnc.pncopen)

Subset on Dimension
~~~~~~~~~~~~~~~~~~~

:func:`~PseudoNetCDF.PseudoNetCDFFile.sliceDimensions` is one of the most basic jobs of PNC. The sliceDimensions method will
extract just the selected elements on one or more dimensions (e.g., LAY). This
method works with `slice` objects, integers, and arrays.

.. code-block:: python

  lay1f = camxf.sliceDimensions(LAY=0)  

Apply Function on a Dimension
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~PseudoNetCDF.PseudoNetCDFFile.applyAlongDimensions` applies functions
on dimensions and is the second most basic jobs of PNC. Functions can be numpy
array method names (e.g., 'mean', 'std', 'sum') or any function.

.. code-block:: python

  timeavgf = camxf.applyAlongDimensions(TSTEP='mean')  


Make a plot
~~~~~~~~~~~

The :func:`~PseudoNetCDF.PseudoNetCDFFile.plot` interface tries to make some
simple plotting easy. It also allows keyword configuration and object-based
interaction with the created plot by returning the axes that was plotted on.

.. code-block:: python

  ax = lay1f.plot('O3')
  ax.figure.savefig('O3.png')

Derive Variables
~~~~~~~~~~~~~~~~

The :func:`~PseudoNetCDF.PseudoNetCDFFile.eval` is very like `pandas.DataFrame.eval`. It requires the non-derived
variables to be in the file or calculated earlier in a multi-line script.

.. code-block:: python

  noxf = lay1f.eval('NOx = NO + NO2')

Stack Files on Dimensions
~~~~~~~~~~~~~~~~~~~~~~~~~

The :func:`~PseudoNetCDF.PseudoNetCDFFile.stack` method allows files to be
concatenated on a dimension. All other dimensions must be of the same length.

.. code-block:: python

  camxpath1 = 'CAMx.v6.40.midwest.36.12.noMPI.20020603.avrg.grd02'
  camxpath2 = 'CAMx.v6.40.midwest.36.12.noMPI.20020604.avrg.grd02'
  camx1f = pnc.pncopen(camxpath1, format='uamiv')
  camx2f = pnc.pncopen(camxpath2, format='uamiv')
  camxf = camx1f.stack(camx2f, 'TSTEP')

Extract Coordinates
~~~~~~~~~~~~~~~~~~~

The :func:`~PseudoNetCDF.PseudoNetCDFFile.ll2ij` converts longitude and
latitude to indices (0-based) and the 
:func:`~PseudoNetCDF.PseudoNetCDFFile.time2t` function converts datetime
objects to time indices (0-based). This makes it easy to extract time and space
elements from a file. Currently, there is not a comparable vertical method.

.. code-block:: python

  from datetime import datetime
  
  obslon = [-100, -90]
  obslat = [30, 40]
  obsdatstrs = [
    datetime(2002, 6, 3, 8, tzinfo=timezone.utc),
    datetime(2002, 6, 3, 9, tzinfo=timezone.utc),
  ]
  
  i, j = camxf.ll2ij(obslon, obslat)
  k = i * 0
  t = camxf.time2t(obstimes)
  
  atobsf = camxf.sliceDimensions(TSTEP=t, LAY=k, ROW=j, COL=i)


Copy a File to Memory
~~~~~~~~~~~~~~~~~~~~~

Finally, it is often nice to make an in-memory copy of a file. Sometimes it can
be as a template or other times it helps methods work more efficiently.

.. code-block:: python

  inmemf = camxf.copy()
