"""
NOAA/GMD sounding data file produced on 2015/02/02 22:04:08 GMT by SkySonde Processor version 0.2.2.8
               Header lines = 45
               Data columns = 21
              Flight number = BL049
                 Date [GMT] = 03-01-1992
                 Time [GMT] = 18:35:21
                   Location = Boulder, CO
                  Longitude = -105.25000
                   Latitude = 40.00000
       Launch altitude (km) = 1.743
      Surface pressure (mb) = 826.790
  Surface temperature (deg) = 6.870
       Surface humidity (%) = 18.900
         Turn altitude (km) = 30.699
         Turn pressure (mb) = 10.20
          Radiosonde number = 8815949
     Vaisala humicap sensor = A
       Radiation correction = Yes
     Pressure sensor offset = 0.000
                 A/D System = 12 bit Tmax
Radiosonde Total Col. Water = 2.322
            Instrument Type = Ozone Sonde
       Original File Source = unknown
         Ozone sonde number = 5A8087
                   Solution = 1%
Oltmans solution correction = Yes
             Oltmans term A = 1.0000
             Oltmans term B = 0.4000
   Total ozone column (CMR) = 283 (41)
  Total ozone column (SBUV) = 289 (47)
  Total ozone stop pressure = 10.28
  Time (sec) to pump 100 ml = 28.230
    Dry flowrate correction = 2.60
    Background current (uA) = 0.071
               Coefficients = 
       Pump coefficient pc0 = 0.5955
                        pc1 = 0.5125
                        pc2 = -0.2353
                        pc3 = 0.0824
       Pressure data source = Radiosonde Pressure
       Altitude data source = Radiosonde Geopot Altitude
   Geo alt anchoring method = First/Launch Row Set to Launch Altitude

      Time,     Press,      Praw,       Alt,     Tcorr,      Temp,      Traw,     Theta,        RH,     RHraw,     TFp V,     IPW V,    TVaisI, O3 Cell I,      O3 P,     O3 Mr,    T Pump,     O3Bat,    I Pump,     Total Column O3,  Total w/ Extrap O3"""
      

from PseudoNetCDF import PseudoNetCDFFile
import re
import numpy as np
spaces = re.compile(r',\s{0,1000}')
def skysonde1sec(inpath):
    datafile = open(inpath, 'r')
    datalines = datafile.read().split('\n')
    nmeta = int(datalines[1].split(' = ')[1])
    meta = dict([[w.strip() for w in l.split(' = ')] for l in datalines[1:nmeta-2] if l != ''])
    varline, unitline = datalines[nmeta-2:nmeta]
    varnames = [vn.strip() for vn in spaces.split(varline)]
    units = [u.strip()[1:-1].strip() for u in spaces.split(unitline)]
    print(units)
    import pdb; pdb.set_trace()
    data = np.fromstring(', '.join(datalines[nmeta:]), sep = ',').reshape(-1, len(varnames))
    outf = PseudoNetCDFFile()
    outf.createDimension('time', data.shape[0])
    for varname, unit, vals in zip(varnames, units, data.T):
        var = outf.createVariable(varname.replace(' ', '_'), 'f', ('time',))
        var.units = unit
        var.standard_name = varname
        var[:] = np.ma.masked_values(vals[:], 99999)
    
    return outf
    
if __name__ == '__main__':
    f = skysonde1sec('/Users/barronh/Downloads/skysonde1sec.txt')
