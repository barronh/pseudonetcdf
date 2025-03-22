.. CMAQ Odds and Ends

Odds and Ends
-------------

This is a collection of small examples that are useful either by themselves
or as examples that can be adapted and combined to more general problems.
The example data used with be IOAPI-like data, but these examples work with
CMAQ, CAMx, GEOS-Chem, HYSPLIT... any file that PseudoNetCDF can read.

Make Test File
~~~~~~~~~~~~~~
The following examples assume the existence of a file (`f`) made here.
which has one variable (DUMMY). Note that this file also has latitude
and longitude as variables that IOAPI would not see.

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
  gf.SDATE = 2016001
  gf.updatetflag(overwrite=True)
  delattr(gf, 'VAR-LIST')
  gf.updatemeta()
  gf.save('test.nc').close()

Create a Variable Like Another
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often variables share dimensions and data types, which makes duplication useful.
The `copyVariable` method will takes a variable and a new `key` to use for the
duplicate. The duplicated variable has all the same properties and data. You
can exclude data copying with `withdata=False` and overwrite the properties.

The example below takes the DUMMY variable and makes a NEW variable. Because
this is IOAPI-like, I remve the VAR-LIST attribute and have it rebuilt.

.. code-block:: python

  import PseudoNetCDF as pnc
  
  f = pnc.pncopen('test.nc', format='ioapi').copy()
  oldvar = f.variables['DUMMY']
  newvar = f.copyVariable(oldvar, key='NEW', withdata=False)
  
  # Update properties
  newvar.setncatts(dict(long_name='NEW', var_desc='NEW variable', units='ppmV'))
  # Update values
  newvar[:] = 1
  delattr(f, 'VAR-LIST')
  f.updatemeta()

Make a Mask
~~~~~~~~~~~
A common task is to create a variable whose shape depends on location in a
geopolitical boundary. This example creates a file with variables for each
United States state. The values are 1 in the state and 0 outside. The political
boundary definitions come from the gadm.org project.

First, we collect shapes for each state. Then, for each cell, we check the
cell centroid against each shape. For efficiency, the state shapes are
prepared (`prep`) and simplied (`simplify`). This allows for fast checking.
Then, we sort the polygons by distance to preferentially test close states.
When a match is found, the loop is broken and moves on to the next cell.

.. code-block:: python
  
  import PseudoNetCDF as pnc
  import numpy as np
  import shapefile
  from shapely.geometry import shape, Point
  from shapely.prepared import prep
  
  
  f = pnc.pncopen('test.nc', format='ioapi').copy()
  tmpvar = f.variables['DUMMY']
  # if your file has lat/lon variables, read them directly
  # or use the code below to make them
  I, J = np.meshgrid(
      np.arange(f.NCOLS),
      np.arange(f.NROWS),
  )
  lon, lat = f.ij2ll(I, J)
  
  # wget https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_USA_shp.zip
  # unzip gadm36_USA_shp.zip
  gadmf = shapefile.Reader('gadm36_USA_1.shp')
  # define useful field ids
  fieldids = {field[0]: fi for fi, field in enumerate(gadmf.fields[1:])}
  statefield = fieldids['NAME_1']
  varfield = fieldids['HASC_1']
  # Create a dictionary of shapes with names like US_NC
  shapes = {}
  for feat in gadmf:
      shapes[feat.record[varfield].replace('.', '_')] = shape(feat.shape)
  
  # make a copy of shapes for fast contains checking
  prepared = {shapek: prep(shape) for shapek, shape in shapes.items()}
  
  # Create variables to hold state masks
  # initialize with 0
  for maskkey in shapes:
      maskvar = f.copyVariable(tmpvar, key=maskkey, withdata=False)
      maskvar.setncatts(dict(units='1', long_name=maskkey, var_desc=maskkey))
      maskvar[:] = 0
  
  # Loop over row (j) and column (i) and find state that contains
  # cell center
  for j in range(f.NROWS):
      print(end='.', flush=True)
      for i in range(f.NCOLS):
          cellcenter = Point(lon[j, i], lat[j, i])
          checkkeys = sorted(list(shapes), key=lambda k: cellcenter.distance(shapes[k].envelope))
          for maskkey in checkkeys:
              prepd = prepared[maskkey]
              if prepd.contains(cellcenter):
                  maskvar = f.variables[maskkey]
                  maskvar[0, 0, j, i] = 1
                  break

  # Add a synthesized variable  
  f.eval("""
  NOAA_NW = US_ID + US_OR + US_WA
  NOAA_NW.long_name = 'NOAA_NW'
  NOAA_NW.var_desc = 'NOAA Northwest Climate Region: Idaho, Oregon and Washingon'
  """, inplace=True)

  # Save as a mask file
  if 'DUMMY' in f.variables:
      del f.variables['DUMMY']
  delattr(f, 'VAR-LIST')
  f.updatemeta()
  f.SDATE = -635
  f.TSTEP = 0
  f.variables['TFLAG'][:, :, 0] = f.SDATE
  f.save('mask.nc')
