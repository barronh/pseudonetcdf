import numpy as np
import re
import io
from PseudoNetCDF import PseudoNetCDFFile
from PseudoNetCDF.pncwarn import warn
from datetime import datetime

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
        try:
            import pandas as pd
        except:
            raise ImportError('woudcsonde requires pandas; install pandas (e.g., pip install pandas)')
        metalines = myf.readline()
        nmetalines = 2
        newblock = ''
        while True:
            nextline = myf.readline()
            nstrip = nextline.strip()
            nmetalines += 1
            if nstrip == '': break
            elif nstrip.startswith('*'): break
            elif nstrip.startswith('#PROFILE'):
                newblock = nstrip[1:] + '_'
                break
            elif nstrip.startswith('#'):
                nmetalines -= 1
                newblock = nstrip.upper()[1:] + '_'
                nlines, newblock = self._addmetavars(myf, key = newblock)
                nmetalines += nlines
                break
            else: metalines += nextline
        
        #try:
        #    locdata = np.recfromcsv(io.BytesIO((metalines).encode()), case_sensitive = True)
        #except ValueError as e:
        if key == 'TIMESTAMP_': dtype = 'S64'
        else: dtype = None
        try:
            locdata = pd.read_csv(io.BytesIO((metalines).encode()), dtype = dtype)
        except Exception as e:
            warn(key + ': ' + str(e))
            return nmetalines
        if key in ('PLATFORM_', 'LOCATION_'):
            dim1 = 'site'
        else:
            dim1 = 'flight'
        for colkey in locdata.columns:#dtype.names:
            val = locdata[colkey]
            varkey = key + colkey 
            dt = val.dtype.char
            odt = {'O': 'c'}.get(dt, dt)
            dim = {'c': (dim1, 'STRLEN'), 'S': (dim1, 'STRLEN')}.get(odt, (dim1,))
            if val.size > 1:
                if not key[:-1] in self.dimensions:
                    self.createDimension(key[:-1], val.size)
                dim = dim + (key[:-1],)
            var = self.createVariable(varkey, odt, dim)
            var.units = _units.get(colkey, 'unknown')
            var.long_name = colkey
            if odt in ('c', 'S'):
                var[:] = np.array([val]).astype('S64').view('S1')
            else:
                var[:] = val
        return nmetalines, newblock
        
    def __init__(self, path, debug = False, na_values = ['\s+', '*', '99999']):
        """
        path - path or file-like object.
        """
        try:
            import pandas as pd
        except:
            raise ImportError('woudcsonde requires pandas; install pandas (e.g., pip install pandas)')
        if hasattr(path, 'readline'):
            myf = path
        else:
            myf = open(path, 'rU')
        
        meta = ''
        li = 0
        self.createDimension('site', 1)
        self.createDimension('flight', 1)
        self.createDimension('STRLEN', 64)
        for l in myf:
            key = l.strip().upper()
            if key == '#PROFILE':
                break
            elif key.startswith('#'):
                if debug: print(key)
                lines, newblock = self._addmetavars(myf, key[1:] + '_')
                li += lines
                if newblock == 'PROFILE_':
                    li -= 1
                    break
            else:
                meta += l
                li += 1
        
        # after #PROFILE loop is broken
        myf.close()
        readopts = dict(sep = '\s*,', skiprows = li + 1, engine = 'python', na_values = na_values, comment = '*')
        data = pd.read_csv(path, **readopts)
        #indkey = data.columns[0]
        indkey = 'level'
        self.createDimension(indkey, len(data))
        self.metadata = meta
        for key  in data.columns:
            try:
                values = data[key]
                var = self.createVariable(key, 'f', ('flight', indkey,))
                var.units = _units.get(key, 'unknown')
                var.long_name = key
                var[:] = values
            except Exception as e:
                warn(str(e) + '; ' + key + ' will not be written')
        
        date = '-'.join(['%02d' % int(v) for v in self.variables['TIMESTAMP_Date'].view('S64')[0,0].decode().strip().split('-')])
        time = (self.variables['TIMESTAMP_Time'].view('S64')[0,0].decode().strip() + ':00')[:8]
        if time == ':00':
            time = '00:00:00'
        z = self.variables['TIMESTAMP_UTCOffset'].view('S64')[0,0].decode().strip()
        if z.startswith('+'): z = z[:3] + '00'
        elif z.startswith('-'): z = z[:3] + '00'
        else: z = '+' + z[:2] + '00'
        datestr = '{} {}{}'.format(date, time, z)
        rdatestr = '1970-01-01 00:00:00+0000'
        rdate = datetime.strptime(rdatestr, '%Y-%m-%d %H:%M:%S%z')
        outdate = datetime.strptime(datestr, '%Y-%m-%d %H:%M:%S%z')
        dt = (outdate - rdate).total_seconds()
        tvar = self.createVariable('time', 'd', ('flight',))
        tvar.long_name = 'time'
        tvar.units = 'seconds since ' + rdatestr
        tvar[...] = dt
    
    def avgSigma(self, vglvls = None, vgtop = None, hyai = None, hybi = None, inplace = False, copyall = True, levelname = 'LAY'):
        """
        vglvls - sigma coordinate
        vgtop - top in Pa

        hyai - hybrid eta-component hPa
        hybi - hybrid sigma-component
        """
        if inplace:
            outf = self
        else:
            outf = self.copy(props = True, dimensions = True, variables = copyall, data = copyall)
        myp = self.variables['Pressure'][:] # hPa
        myo = self.variables['O3PartialPressure'][:] * 1e1 # mPa -> cPa
        psfc = myp[0, 0]
        if not vglvls is None:
            outpress_edges = vglvls * (psfc - vgtop / 100) + vgtop / 100
            nlvls = vglvls.size - 1
        elif not hyai is None:
            outpress_edges = hyai + hybi * psfc
            nlvls = hyai.size - 1
        outf.createDimension(levelname, nlvls)
        abovebottom = myp[0, :, None] <= outpress_edges[None, :-1]
        belowtop = myp[0, :, None] > outpress_edges[None, 1:]
        invglvl = abovebottom & belowtop
        out = np.ma.masked_all((1, nlvls), dtype = 'd')
        for li in range(nlvls):
            idx = invglvl[:, li]
            if not idx.any(): continue
            lps = myp[:, idx]
            lops = np.ma.masked_invalid(myo[:, idx])
            lopsmask = np.ma.getmaskarray(lops)
            lp = np.ma.masked_where(lopsmask, lps).sum(-1)
            lop = lops.sum(-1)
            out[:, li] = lop / lp
        outv = outf.createVariable('O3', 'f', ('flight', levelname), fill_value = -999)
        outv.units = 'ppm'
        outv.missing_value = -999
        outv.description = 'Pressure-weighted average vmr in ppm'
        outv[:] = out[:]
        return outf

if __name__ == "__main__":
    f = woudcsonde('20150813.ECC.1Z.1Z27925.JMA.csv')
