__all__ = ['arlpdump']
from PseudoNetCDF import PseudoNetCDFFile
import numpy as np
class arlpdump(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, path):
        try:
            f = arlpdump(path)
            return True
        except:
            return False
    def __init__(self, path):
        self._path = path
        self._file = open(path, 'rb')
        self._read()
        mass = self._data['data']['f1']
        loc = self._data['data']['f4']
        meta = self._data['data']['f7']
        self.createDimension('PARTICLES', self.NPARTICLES)
        self.createDimension('POLLUTANTS', self.NPOLLUTANTS)
        pvar = self.createVariable('PARTICLE_MASS', 'f', ('PARTICLES', 'POLLUTANTS'), values = mass)
        pvar.units = 'mass'
        pvar.long_name = 'PARTICLE_MASS'
        var = self.createVariable('latitude', 'f', ('PARTICLES',), values = loc[:,0])
        var.units = 'degrees_north'
        var.long_name = 'latitude'
        var = self.createVariable('longitude', 'f', ('PARTICLES',), values = loc[:,1])
        var.units = 'degrees_east'
        var.long_name = 'latitude'
        var = self.createVariable('height', 'f', ('PARTICLES',), values = loc[:,2])
        var.units = 'meters'
        var.long_name = 'height'
        var = self.createVariable('sigma_u', 'f', ('PARTICLES',), values = loc[:,3])
        var.units = 'sigma'
        var.long_name = 'SIGMA-U'
        var = self.createVariable('sigma_v', 'f', ('PARTICLES',), values = loc[:,4])
        var.units = 'sigma'
        var.long_name = 'SIGMA-V'
        var = self.createVariable('sigma_x', 'f', ('PARTICLES',), values = loc[:,5])
        var.units = 'sigma'
        var.long_name = 'SIGMA-X'
        var = self.createVariable('AGE', 'i', ('PARTICLES',), values = meta[:,0])
        var.units = 'seconds?'
        var.long_name = 'AGE'
        var = self.createVariable('DISTRIBUTION', 'i', ('PARTICLES',), values = meta[:,1])
        var.units = 'unknown'
        var.long_name = 'DISTRIBUTION'
        var = self.createVariable('POLLUTANT', 'i', ('PARTICLES',), values = meta[:,2])
        var.units = 'unknown'
        var.long_name = 'POLLUTANT'
        var = self.createVariable('METEO_GRID', 'i', ('PARTICLES',), values = meta[:,3])
        var.units = 'unknown'
        var.long_name = 'METEO-GRID'
        var = self.createVariable('SORT_INDEX', 'i', ('PARTICLES',), values = meta[:,4])
        var.units = 'none'
        var.long_name = 'SORT-INDEX'
    
    def _read(self):
        """
        INT*4 Number of particles
        INT*4 Number of pollutants
        INT*4 Time of particle dump (YEAR, MONTH, DAY, HOUR, MINUTES)
        """
        rec1 = np.fromfile(self._path, dtype = '>i,>7i,>i', count = 1)
        self._file.seek(0, 0)
        hdr = rec1['f1'][0]
        self.NPARTICLES = hdr[0]
        self.NPOLLUTANTS = hdr[1]
        self.YEAR = hdr[2]
        self.MONTH = hdr[3]
        self.DAY = hdr[4]
        self.HOUR = hdr[5]
        self.MINUTES = hdr[6]
        assert(rec1['f0'] == rec1['f2'])
        self._fmt = np.dtype([('header', '>i,>7i,>i'), ('data', '>i,>({},)f,>i,>i,>6f,>i,>i,>5i,>i'.format(self.NPOLLUTANTS), (self.NPARTICLES,))]);
        self._data = np.fromfile(self._path, dtype = self._fmt, count = 1)[0]
        db = self._data['data']
        assert((db['f0'] == db['f2']).all())
        assert((db['f3'] == db['f5']).all())
        assert((db['f6'] == db['f8']).all())
    
if __name__ == '__main__':
    from PseudoNetCDF import pncopen
    f = pncopen('PARDUMP_143', format = 'arlpdump')
