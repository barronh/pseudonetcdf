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

  # Either open GRIDDESC from disk or create one in memory
  gdpath = '/path/to/GRIDDESC'
  gdpath = """' '
  'LamCon_40N_97W'\n 2 33.000 45.000 -97.000 -97.000 40.000
  ' '
  '12US1'
  'LamCon_40N_97W' -2556000.0 -1728000.0 12000.0 12000.0 459 299 1
  ' '
  """
  gf = pnc.pncopen(gdpath, format='griddesc', FTYPE=2, SDATE=1970001, TSTEP=0)
 
  # Assumes you have a CMAQ 3D CONC file on same vertical grid 
  concf = pnc.pncopen('CONC.nc', format='ioapi')
  i, j = concf.ll2ij(
      gf.variables['longitude'][:],
      gf.variables['latitude'][:]
  )
  bconf = concf.sliceDimensions(COL=i, ROW=j, newdims=('PERIM',))
  bconf.FTYPE = 2
  bconf.save('bcon.nc', format='NETCDF3_CLASSIC')

