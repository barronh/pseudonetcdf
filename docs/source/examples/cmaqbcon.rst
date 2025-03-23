.. CMAQ Simple Boundary Conditions

CMAQ Simple Boundary Conditions
-------------------------------

Creating boundary conditions can be complex and is solved more completely
using github.com/barronh/aqmbc or the BCON process in github.com/USEPA/CMAQ.
This example extracts perimeter locations and applies no other processing.
It is meant to be an example of capability, not a replacement for those
tools.

The example assumes the existence of an IOAPI file that is bigger than the
domain for which boundaries are being extracted. `CONC.nc` is from a H-CMAQ
simulation, but could have come from a coarser domain.

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
  gf = pnc.pncopen('GRIDDESC', format='griddesc', FTYPE=2)
  concf = pnc.pncopen('CONC.nc', format='ioapi')
  i, j = concf.ll2ij(
      gf.variables['longitude'][:],
      gf.variables['latitude'][:]
  )
  bconf = concf.sliceDimensions(COL=i, ROW=j, newdims=('PERIM',))
  bconf.FTYPE = 2
  bconf.save('bcon.nc', format='NETCDF3_CLASSIC')

