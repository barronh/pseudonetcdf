__all__ = ['l100']
from datetime import datetime
from PseudoNetCDF import PseudoNetCDFFile
import numpy as np

_edges = [0, 4, 12, 20, 28, 35, 43, 50, 57, 64, 72, 79, 88, 95]
_starts = _edges[:-1]
_ends = _edges[1:]
_ses = [se for se in zip(_starts, _ends)]
_orignames = 'Level Press Alt Pottp Temp FtempV Hum Ozone Ozone Ozone Ptemp O3_DN O3_Res'.split()
_names     = 'Level Press Alt Pottp Temp FtempV Hum Ozone_mPa Ozone_ppmv Ozone_atmcm Ptemp O3_DN O3_Res'.split()
_units     = 'Num hPa km K C C % mPa ppmv atmcm C 10^11/cc DU'.split()


class l100(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, path):
        """
        Must have correct names on line 28
        """
        try:
            lines = cls._getmeta(path)
            mynames = lines[-2].split()
            for chk, new in zip(_orignames[:8], mynames):
                if chk != new: return False
            else:
                return True
        except:
            return False
    
    @classmethod
    def _getmeta(cls, path):
        
        # get metadata lines
        infile = open(path)
        lines = []
        for li in range(100):
            line = infile.readline()
            lines.append(line)
            if line.startswith('Level'):
                lines.append(infile.readline())
                break
        else:
            lines = lines[:29]
        return lines
        
    def __init__(self, path):
        try:
            import pandas as pd
        except:
            raise ImportError('l100 sonde files requires pandas; install pandas (e.g., pip install pandas)')
        self._path = path
        metalines = l100._getmeta(path)
        datastartline = len(metalines)
        for line in metalines:
            if ':' in line:
                parts = line.split(':')
                key = parts[0].strip().replace(' ', '_')
                val = ':'.join(parts[1:]).strip()
                setattr(self, key, val)
        
        dataset = pd.read_fwf(path, colspecs = _ses, skiprows = datastartline, names = _names, na_values = ['99.90', '99.99', '999', '99.999', '999.999', '99.9990', '9999'])
        self.createDimension('site', 1)
        self.createDimension('level', dataset.shape[0])
        for key, unit in zip(_names, _units):
            vals = dataset[key]
            outvar = self.createVariable(key, vals.dtype, ('site', 'level'))
            outvar.long_name = key
            outvar.units = unit
            outvar[0, :] = vals
        
        ldate = datetime.strptime(self.Launch_Date + ' ' + self.Launch_Time, '%d %B %Y %H:%M:%S %Z')
        rdate = datetime(1970,1,1)
        seconds = (ldate - rdate).total_seconds()
        lonvar = self.createVariable('longitude', 'd', ('site',))
        lonvar.long_name = 'longitude'
        lonvar.units = 'degrees'
        lonvar[:] = eval(self.Longitude)
        latvar = self.createVariable('latitude', 'd', ('site',))
        latvar.long_name = 'latitude'
        latvar.units = 'degrees'
        latvar[:] = eval(self.Latitude)
        tvar = self.createVariable('time', 'd', ('site',))
        tvar.units = 'seconds since 1970-01-01 00:00:00+0000'
        tvar[:] = seconds
    
        self.description = ''.join(metalines)
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
        myp = self.variables['Press'][:] # hPa
        myo = self.variables['Ozone_mPa'][:] * 1e1 # mPa -> cPa
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
        outv = outf.createVariable('O3', 'f', ('site', levelname), missing_value = -999)
        outv.units = 'ppm'
        outv.description = 'Pressure-weighted average vmr in ppm'
        outv[:] = np.ma.filled(out[:], np.nan)
        return outf

