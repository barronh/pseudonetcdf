import pandas as pd
import numpy as np
from PseudoNetCDF import PseudoNetCDFFile

_units = dict(
          Pressure = 'hPa', O3PartialPressure = 'mPa',
          Temperature = 'C', WindSpeed = 'm/s', WindDirection = 'degrees from north', RelativeHumidity = '%',
          LevelCode = '?', Duration = 's', GPHeight = 'm', SampleTemperature = 'C',
          Latitude = 'degrees_north', Longitude = 'degrees_east')
        

class woudcsonde(PseudoNetCDFFile):
    """wouldcsonde provides a NetCDF-like interface to csv files from the 
    ozonesonde data distributed by the Would Ozone and Ultraviolet
    Radiation Data Center (wouldc)
    
    http://woudc.org/data/products/
    """
    
    def _addmetavars(self, myf, key = ''):
        """Adds metavariables for next two lines"""
        metalines = myf.readline()
        nmetalines = 2
        while True:
            nextline = myf.readline()
            nmetalines += 1
            if nextline.strip() == '': break
            elif nextline.strip().startswith('*'): break
            else: metalines += nextline
        import io
        try:
            locdata = np.recfromcsv(io.BytesIO((metalines).encode()), case_sensitive = True)
        except:
            return nmetalines
        self.createDimension('site', 1)
        self.createDimension('STRLEN', 64)
        for colkey in locdata.dtype.names:
            val = locdata[colkey]
            varkey = key + colkey 
            dt = val.dtype.char
            dim = {'c': ('site', 'STRLEN'), 'S': ('site', 'STRLEN')}.get(dt, ('site',))
            if val.size > 1:
                if not key in self.dimensions:
                    self.createDimension(key, val.size)
                dim = dim + (key,)
            var = self.createVariable(varkey, dt, dim)
            var.units = _units.get(colkey, 'unknown')
            var.long_name = colkey
            if dt in ('c', 'S'):
                var[:] = np.array([val]).astype('S64').view('S1')
            else:
                var[:] = val
        return nmetalines
        
    def __init__(self, path, debug = False):
        """
        path - path or file-like object.
        """
        if hasattr(path, 'readline'):
            myf = path
        else:
            myf = open(path, 'rU')
        
        meta = ''
        li = 0
        for l in myf:
            key = l.strip().upper()
            if key == '#PROFILE':
                break
            elif key.startswith('#'):
                if debug: print(key)
                lines = self._addmetavars(myf, key[1:] + '_')
                li += lines
            else:
                meta += l
                li += 1
        
        # after #PROFILE loop is broken
        myf.close()
        data = pd.read_csv(path, sep = ',', skiprows = li + 1)
        indkey = data.columns[0]
        self.createDimension(indkey, len(data))
        self.metadata = meta
        for key  in data.columns:
            values = data[key]
            var = self.createVariable(key, 'f', ('site', indkey,))
            var.units = _units.get(key, 'unknown')
            var.long_name = key
            var[:] = values
    
    def avgSigma(self, vglvls, vgtop, inplace = False, copyall = True):
        """
        vglvls - sigma coordinate
        vgtop - top in Pa
        """
        if inplace:
            outf = self
        else:
            outf = self.copy(props = True, dimensions = True, variables = copyall, data = copyall)
        nlvls = vglvls.size - 1
        outf.createDimension('LAY', nlvls)
        myp = self.variables['Pressure'][:] * 100 # hPa -> Pa
        myo = self.variables['O3PartialPressure'][:] * 1e3 # mPa -> uPa
        psfc = myp[0, 0]
        ptop = vgtop
        mysigma =  (myp - ptop) / (psfc - ptop)
        abovebottom = mysigma[0, :, None] <= vglvls[None, :-1]
        belowtop = mysigma[0, :, None] > vglvls[None, 1:]
        invglvl = abovebottom & belowtop
        out = np.ma.masked_all((1, nlvls), dtype = 'd')
        for li in range(nlvls):
            idx = invglvl[:, li]
            if not idx.any(): continue
            lp = myp[:, idx].sum(-1)
            lop = myo[:, idx].sum(-1)
            out[:, li] = lop / lp
        outv = outf.createVariable('O3', 'f', ('site', 'LAY'))
        outv.units = 'ppm'
        outv.description = 'Pressure-weighted average vmr in ppm'
        outv[:] = out[:]
        return outf
        
if __name__ == "__main__":
    f = woudcsonde('20150813.ECC.1Z.1Z27925.JMA.csv')
