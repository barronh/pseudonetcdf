"""
.. CMAQ Spatial Subset

CMAQ Spatial Subset
-------------------

This example shows how to spatially subset a file in two ways. The first
subset windows a file and retains the source file resolution. The second
subset resamples the coarse data using cell centroid to extract concentrations
from the source file and put it in a finer resolution file.

Both examples use he extent of the 12US1 file defined in GRIDDESC to create
a bounding box and apply it to the 108NHEMI2 domain data (`CONC.nc`).

CMAQ Window
~~~~~~~~~~~

This first code block takes a polar stereographic file and extracts just
the cells within the 12US1 coninental domain. The meta-data is updated
automatically, so this is pretty easy.

.. code-block:: python

  import PseudoNetCDF as pnc
  
  # create a GRIDDESC
  with open('GRIDDESC', 'w') as gf:
      gf.write(
         "' '\n'LamCon_40N_97W'\n 2 33.000 45.000 -97.000 -97.000 40.000\n" +
         "' '\n'12US1'\n'LamCon_40N_97W' " +
         "-2556000.0 -1728000.0 12000.0 12000.0 459 299 1\n' '"
      )

  # Create a Template
  gf = pnc.pncopen('GRIDDESC', format='griddesc')
  concf = pnc.pncopen('CONC.nc', format='ioapi').copy()
  
  i, j = concf.ll2ij(
      gf.variables['longitude'][:],
      gf.variables['latitude'][:]
  )
  islice = slice(i.min(), i.max() + 1)
  jslice = slice(j.min(), j.max() + 1)
  wndwf = concf.sliceDimensions(COL=islice, ROW=jslice)
  wndwf.save('wndw.nc', format='NETCDF3_CLASSIC')


CMAQ Resampling
~~~~~~~~~~~~~~~

This second code block is a continuation of the one above. Instead of
windowing and preserving source resolution, we are resampling to produce
a file of the target resolution. The extraction method is complex and so
much of the code is spent making sure the meta-data is kept up to date.

.. code-block:: python

  # Get values for each i/j
  tmpf = concf.sliceDimensions(LAY=0).sliceDimensions(COL=i, ROW=j, newdims=('ROW', 'COL'))

  startdate = concf.getTimes()[0]
  # Create an output file for resampled data
  gf.updatetflag(startdate=startdate, overwrite=True)
  gf.VGLVLS = tmpf.VGLVLS
  nsteps = len(tmpf.dimensions['TSTEP'])
  resamplef = gf.subsetVariables(['DUMMY']).sliceDimensions(
      TSTEP=[0]*nsteps
  )
  # Set times to be sequential
  resamplef.updatetflag(startdate=startdate, tstep=concf.TSTEP, overwrite=True)
  # Copy in the data
  for key, var in tmpf.variables.items():
      if key == 'TFLAG': continue
      resamplef.copyVariable(var, key=key)
  
  del resamplef.variables['DUMMY']
  delattr(resamplef, 'VAR-LIST')
  resamplef.updatemeta()
  # Make sure TFLAG matches 
  resamplef.updatetflag(tstep=concf.TSTEP, overwrite=True)
  # Save file
  resamplef.save('resample.nc')
