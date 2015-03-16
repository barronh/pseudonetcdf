# simple\_xy #

Following the NetCDF tutorial, we'll create a simple file with 2 dimensions (x,y) and one variable.  Comments begin with # and are unnecessary.

```
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from PseudoNetCDF.pncdump import pncdump
from numpy import arange

# create a file like object
simple_xy_file = PseudoNetCDFFile()

# create the x dimension
simple_xy_file.createDimension('x', 6)

# create the x dimension
simple_xy_file.createDimension('y', 12)

# create the data varaible
data = simple_xy_file.createVariable('data', 'i', ('x', 'y'))

# create a numpy array with the right data
data_vals = arange(72).reshape(6, 12)

# load the data
data[:] = data_vals

# dump the file using the PseudoNetCDF
# analog to ncdump
pncdump(simple_xy_file, name = 'simple_xy')
```

You should see a file almost exactly like the one in the [NetCDF tutorial](http://www.unidata.ucar.edu/software/netcdf/docs_beta/netcdf-tutorial.html#simple_005fxy).  You can even pipe the output of pncdump to ncgen and it will make a real netcdf file.