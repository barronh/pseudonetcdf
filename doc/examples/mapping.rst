.. Mapping

Mapping Explained
-----------------

The `PseudoNetCDF.plot` function makes maps when `plottype` is `longitude-latitude`, which is the default (in new versions, the default is the last two dimensions; often COL-ROW, lon-lat, x-y). When `plot` triggers a map, it can be thought of in four parts:

1. Get a map
2. Draw features
3. Add a pseudocolor (aka tile) feature from the variables
4. Add a colorbar

To illustrate, how it works the three following code segments perform the same task using successively less code. In all three examples, the ACONC.nc file is an IOAPI output from CMAQ that contains the NO2 variable.

Mapping with Basemap
~~~~~~~~~~~~~~~~~~~~

The first segment uses `netCDF4` to open the file, `numpy` to construct coordinates, `basemap` to create a map, and `matplotlib` to do the plotting. 

.. code-block:: python

  import matplotlib.pyplot as plt
  import numpy as np
  import matplotlib.colors as mc
  from mpl_toolkits.basemap import Basemap
  from netCDF4 import Dataset
  
  path = 'ACONC.nc'
  f = Dataset(path)
  var = f.variables['NO2']
  
  bmap = Basemap(
      projection='lcc',
      lon_0=f.P_GAM, lat_0=f.YCENT, lat_1=f.P_ALP, lat_2=f.P_BET,
      llcrnrx=f.XORIG, llcrnry=f.YORIG,
      urcrnrx=f.XORIG + f.NCOLS * f.XCELL,
      urcrnry=f.YORIG + f.NROWS * f.YCELL
  )
  
  x = np.arange(f.NCOLS + 1) * f.XCELL
  y = np.arange(f.NROWS + 1) * f.YCELL
  
  fig, ax = plt.subplots(1, 1)
  bmap.drawcoastlines(ax=ax)
  bmap.drawcountries(ax=ax)
  lay1avgtime = var[:, 0].mean(0)
  
  plt.pcolormesh(x, y, lay1avgtime, norm=mc.LogNorm())
  plt.colorbar(label=varkey + ' ' + var.units.strip())

Mapping Partial PseudoNetCDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second code segment uses `PseudoNetCDF` to open and get the map, but is otherwise similar to the first.

.. code-block:: python

  import matplotlib.pyplot as plt
  import numpy as np
  import matplotlib.colors as mc
  import PseudoNetCDF as pnc
  
  path = 'ACONC.nc'
  f = pnc.pncopen(path, format='ioapi')
  var = f.variables['NO2']
  
  bmap = f.getMap()
  
  x = np.arange(f.NCOLS + 1) * f.XCELL
  y = np.arange(f.NROWS + 1) * f.YCELL
  
  fig, ax = plt.subplots(1, 1)
  bmap.drawcoastlines(ax=ax)
  bmap.drawcountries(ax=ax)
  lay1avgtime = var[:, 0].mean(0)
  
  plt.pcolormesh(x, y, lay1avgtime, norm=mc.LogNorm())
  plt.colorbar(label=varkey + ' ' + var.units.strip())


Mapping with PseudoNetCDF Only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This third version relies on `PseudoNetCDF` to perform all steps.

.. code-block:: python

  import matplotlib.colors as mc
  import PseudoNetCDF as pnc
  
  path = 'ACONC.nc'
  f = pnc.pncopen(path, format='ioapi')
  ax = f.plot('NO2', plot_kw=dict(norm=mc.LogNorm()))


Additional Options for PseudoNetCDF.plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The `PseudoNetCDF.plot` function takes `map_kw`, `plot_kw`, and `cbar_kw` options that let you further configure the plot.

* `map_kw` takes a dictionary that is *mostly* passed through to the `PseudoNetCDF.getMap` function. In addition, it takes a boolean for for each of the following `coastlines`, `countries`, `states`, and `counties`. Setting these to `True` will add that feature to the map. At this time `counties` is automatically set to a low z-order and may not be visible.
* `plot_kw` takes a dictionary that is passed through as keywords to the `matplotlib.Axes.pcolormesh` method.
* `cbar_kw` takes a dictionary that passed through as keywords to the  `matplotlib.Figure.colorbar` method.

