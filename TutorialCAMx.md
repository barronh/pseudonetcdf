# CAMx File Reader Tutorial #

from PseudoNetCDF.camxfiles.Memmaps import 

Gets you input and output (memmory map) readers for CAMx file formats.  Different CAMx files have different levels of meta data; the metadata ranges from completely absent  to nearly complete (i.e. just missing projection).  The  output formats have pretty good meta data.  Some of
the input formats (e.g. vertical\_diffusivity or kv and height\_pressure) do not have sufficient metadata to define the horizontal dimensions (i.e. NROWS x  NCOLS).  Where there is insufficient meta data, there are generally optional rows and cols integer arguments.  If they are not specified, rows will be set to rows **cols and cols is set to 1.  In some cases, the rows and cols arguments are required.  The various readers are listed below with a quick description of meta  data.**

output classes:
  * uamiv: UAM-IV formatted gridded files (initial concentration, gridded emission, etc)
  * ipr: Integrated Process Rate files
  * irr: Integrated Reaction Rate files

inputs classes:
  * humidity: requires rows and cols
  * landuse: requires rows and cols
  * wind: requires rows and cols
  * vertical\_diffusivity: requires rows and cols
  * height\_pressure: optional, but important rows and cols
  * temperature: optional, but important rows and cols
  * cloud\_rain: no required arguments
  * point\_source: no required arguments
  * uamiv: no required arguments; nearly complete metadata; reads


For any of the readers above, you can initialize a file by typing the class name followed by the path and, when appropriate, rows and cols.
For more help on any, type help(classname).

The example below opens one of each input file and an output average conc file for the 2005 May TCEQ simulation using the tceq2005ep0 meteorology, the reg8 gridded emissions, and the pscfv2 point sources.  After opening each file, the script performs example calculations.

```
# import libraries
from PseudoNetCDF.camxfiles.Memmaps import *
from numpy import *

# set file paths
windpath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_wind.20050519.hgbpa_04km.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.v45'
zppath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_zp.20050519.hgbpa_04km.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.v45'
kvpath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_kv.20050519.hgbpa_04km.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.TKE.v45'
lupath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_landuse.hgbpa_04km.utcsr'
temppath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_temp.20050519.hgbpa_04km.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.v45'
crpath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_cr.20050519.hgbpa_04km.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.v45'
humpath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0/grid_04k/camx_hum.20050519.hgbpa_04km.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.v45'
ptpath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0_reg8_pscfv2/grid_04k/camx_cb05_ei_el.20050519.hgb8h2.bc05may.reg7a_pscfv2'
gepath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/inputs/CAMx/tceq2005ep0_reg8_pscfv2/grid_04k/camx_cb05_ei_lo.20050519.hgb8h2.bc05may.reg8.hgbpa_04km'
aconcpath = '/Volumes/ServerRAID/simulations/2005_05_hg_TCEQ/predicted/CAMx/tceq2005ep0_reg8_pscfv2/grid_04k/camx451_cb05.20050519.hgb8h2.bc05may.reg8_pscfv2.2005ep0_eta_dbemis_fddats_uhsst_utcsrlulc.avrg.hgbpa_04km'

# Open files
windfile = wind(windpath, rows = 65, cols = 83)
zpfile = height_pressure(zppath, rows = 65, cols = 83)
kvfile = vertical_diffusivity(kvpath, rows = 65, cols = 83)
lufile = landuse(lupath, rows = 65, cols = 83)
tempfile = temperature(temppath, rows = 65, cols = 83)
crfile = cloud_rain(crpath, rows = 65, cols = 83)
humfile = humidity(humpath, rows = 65, cols = 83)
ptfile = point_source(ptpath)
gefile = uamiv(gepath)
aconcfile = uamiv(aconcpath)

# Perform demonstration calculations
print "Gridded emission file variable keys", gefile.variables.keys()
print "Gridded NO emission rate unit", gefile.variables['NO'].units
print "Gridded emission NOx average over space in %s" % gefile.variables['NO'].units, eval('NO[:] + NO2[:]', globals(), gefile.variables).sum(3).sum(2).sum(1)
print "Maximum pressure in %s" % zpfile.variables['PRES'].units, zpfile.variables['PRES'][:].max()
print "Average layer 1, noon temperature in %s" % tempfile.variables['AIRTEMP'].units, tempfile.variables['AIRTEMP'][12, 0].mean()
print "Maximum cloud water content in %s" % crfile.variables['CLOUD'].units, crfile.variables['CLOUD'][:].mean()
print "Median wind speed in %s" % windfile.variables['U'].units, median(eval('sqrt(U**2 + V**2)', globals(), windfile.variables))
print "Ozone Changes larger than 20ppb/hr", eval('diff(O3, axis = 0) > 0.02', globals(), aconcfile.variables).sum()

R = 831.4472 #  m**3 hPa K**-1 mol**-1
Na = 6.02214179e23 # avagadro's number molecules/mol
dZ = eval('concatenate([HGHT[:,[0]], diff(HGHT[:], axis = 1)], axis = 1)', globals(), zpfile.variables) # m
vol = dZ * gefile.XCELL * gefile.YCELL # m**3
temp = tempfile.variables['AIRTEMP'][:] # K
pres = zpfile.variables['PRES'][:] # hPa
print "Min air number density (molecules/m**3)", (pres * vol / temp / R * Na).min()
## --- END EXAMPLE
```