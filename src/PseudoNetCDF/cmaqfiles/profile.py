__all__ = ['bcon_profile', 'icon_profile']
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariables
from matplotlib.mlab import csv2rec
from StringIO import StringIO
from datetime import datetime
import numpy as np
class icon_profile(PseudoNetCDFFile):
    def __init__(self, path):
        lines = open(path).read().split('\n')
        header = lines[3].split()
        nlay, nspc = map(int, header[:2])
        sigmas = map(float, header[2:])
        nsigmas = len(sigmas)
        try:
            dateo = datetime.strptime(lines[4].strip(), '%Y-%m-%d')
        except:
            date, time = map(int, lines[4].split())
        starts =  [5]
        ends = [s + nspc for s in starts]
        keys = ['all']
        fieldnames = ('name',) + tuple(['s%f' % i for i in sigmas])
        data = dict([(k, csv2rec(StringIO('\n'.join(lines[s:e])), delimiter = ' ', names = fieldnames, converterd = dict(names = lambda x: str(x).strip()))) for k, s, e in zip(keys, starts, ends)])
        profile_spcs = np.char.strip(data[keys[0]]['name'])
        data_type = data[keys[0]].dtype
        data_shape =  data[keys[0]].shape
        self.createDimension('sigma', nsigmas)
        self.createDimension('sigma-mid', nsigmas - 1)
        self.createDimension('south_east_north_west', 4)
        self.createVariable('sigma', 'f', ('sigma',), values = np.array(sigmas), units = 'sigma')
        self.createVariable('sigma-mid', 'f', ('sigma-mid',), values = np.array(sigmas).repeat(2, 0)[1:-1].reshape(-1, 2).mean(1), units = 'sigma')
        self.VGLVLS = self.variables['sigma']
        self.VGTOP = 5000
        ks = keys[1:]
        for k in ks:
            try:
                assert((np.char.strip(data[k]['name']) == profile_spcs).all())
                assert(data[k].dtype == data_type)
                assert(data[k].dtype == data_type)
            except AssertionError:
                raise IOError('File is corrupt or inconsistent')
        for a in data['all']:
            self.createVariable(a[0].strip(), 'f', ('sigma-mid', 'south_east_north_west'), units = "None", values = np.array(map(lambda x: tuple(x)[1:], [a])).T, long_name = a[0].ljust(16), var_desc = a[0].ljust(16))

class bcon_profile(PseudoNetCDFFile):
    def __init__(self, path):
        lines = open(path).read().split('\n')
        header = lines[3].split()
        nlay, nspc = map(int, header[:2])
        sigmas = map(float, header[2:])
        nsigmas = len(sigmas)
        try:
            dateo = datetime.strptime(lines[4].strip(), '%Y-%m-%d')
        except:
            date, time = map(int, lines[4].split())
        starts =  [5 + i + i * nspc for i in range(4)]
        ends = [s + 1 + nspc for s in starts]
        keys = [lines[s].strip().lower() for s in starts]
        fieldnames = ('name',) + tuple(['s%f' % i for i in sigmas])
        data = dict([(k, csv2rec(StringIO('\n'.join(lines[s+1:e])), delimiter = ' ', names = fieldnames, converterd = dict(names = lambda x: str(x).strip()))) for k, s, e in zip(keys, starts, ends)])
        profile_spcs = np.char.strip(data[keys[0]]['name'])
        data_type = data[keys[0]].dtype
        data_shape =  data[keys[0]].shape
        self.createDimension('sigma', nsigmas)
        self.createDimension('sigma-mid', nsigmas - 1)
        self.createDimension('south_east_north_west', 4)
        self.createVariable('sigma', 'f', ('sigma',), values = np.array(sigmas), units = 'sigma')
        self.createVariable('sigma-mid', 'f', ('sigma-mid',), values = np.array(sigmas).repeat(2, 0)[1:-1].reshape(-1, 2).mean(1), units = 'sigma')
        ks = keys[1:]
        self.VGLVLS = self.variables['sigma']
        self.VGTOP = 5000
        for k in ks:
            try:
                assert((np.char.strip(data[k]['name']) == profile_spcs).all())
                assert(data[k].dtype == data_type)
                assert(data[k].dtype == data_type)
            except AssertionError:
                raise IOError('File is corrupt or inconsistent')
        for w, s, e, n in zip(data['west'], data['south'], data['east'], data['north']):
            assert(w[0] == s[0] and e[0] == n[0] and n[0] == s[0])
            self.createVariable(w[0].strip(), 'f', ('sigma-mid', 'south_east_north_west'), units = "None", values = np.array(map(lambda x: tuple(x)[1:], [s, e, n, w])).T, long_name = w[0].ljust(16), var_desc = w[0].ljust(16))

if __name__ == '__main__':
    po = profile('testdata/profile.dat')
    from PseudoNetCDF.pncdump import pncdump
    pncdump(po)
    print po.variables['ATOL1J'].shape
    print po.variables['ATOL1J'].dimensions