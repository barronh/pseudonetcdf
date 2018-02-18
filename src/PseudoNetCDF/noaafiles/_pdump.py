__all__ = ['arlpardump']
from PseudoNetCDF import PseudoNetCDFFile
import numpy as np
from datetime import datetime
class arlpardump(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, path):
        try:
            f = arlpardump(path)
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
        self.createDimension('particles', self.NPARTICLES)
        self.createDimension('pollutants', self.NPOLLUTANTS)
        startdstr = '{:02d}{:02d}{:02d}{:02d}{:02d}+0000'.format(self.YEAR, self.MONTH, self.DAY, self.HOUR, self.MINUTES)
        startd = datetime.strptime(startdstr, '%y%m%d%H%M%z')
        pvar = self.createVariable('particle_mass', 'f', ('particles', 'pollutants'), values = mass)
        pvar.units = 'arbitrary'
        pvar.long_name = 'PARTICLE_MASS'
        var = self.createVariable('latitude', 'f', ('particles',), values = loc[:,0])
        var.units = 'degrees_north'
        var.long_name = 'latitude'
        var = self.createVariable('longitude', 'f', ('particles',), values = loc[:,1])
        var.units = 'degrees_east'
        var.long_name = 'latitude'
        var = self.createVariable('height', 'f', ('particles',), values = loc[:,2])
        var.units = 'meters'
        var.long_name = 'HEIGHT'
        var = self.createVariable('sigma_u', 'f', ('particles',), values = loc[:,3])
        var.units = 'sigma'
        var.long_name = 'SIGMA-U'
        var = self.createVariable('sigma_v', 'f', ('particles',), values = loc[:,4])
        var.units = 'sigma'
        var.long_name = 'SIGMA-V'
        var = self.createVariable('sigma_x', 'f', ('particles',), values = loc[:,5])
        var.units = 'sigma'
        var.long_name = 'SIGMA-X'
        var = self.createVariable('age', 'i', ('particles',), values = meta[:,0])
        var.units = 'minutes since {}'.format(startd.strftime('%Y-%m-%d %H:%M:%S%z'))
        var.long_name = 'AGE'
        var = self.createVariable('distribution', 'i', ('particles',), values = meta[:,1])
        var.units = '---'
        var.long_name = 'DISTRIBUTION'
        var = self.createVariable('pollutant', 'i', ('particles',), values = meta[:,2])
        var.units = '---'
        var.long_name = 'POLLUTANT'
        var = self.createVariable('meteo_grid', 'i', ('particles',), values = meta[:,3])
        var.units = '---'
        var.long_name = 'METEO-GRID'
        var = self.createVariable('sort_index', 'i', ('particles',), values = meta[:,4])
        var.units = '---'
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
        data = np.fromfile(self._path, dtype = self._fmt)
        assert(data.shape[0] == 1)
        self._data = data[0]
        db = self._data['data']
        assert((db['f0'] == db['f2']).all())
        assert((db['f3'] == db['f5']).all())
        assert((db['f6'] == db['f8']).all())
    
if __name__ == '__main__':
    from PseudoNetCDF import pncopen
    f = pncopen('PARDUMP_143', format = 'arlpardump')
