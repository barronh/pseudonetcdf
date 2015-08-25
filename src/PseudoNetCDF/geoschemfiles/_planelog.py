import numpy as np
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable
from glob import glob
from csv import DictReader, Sniffer, QUOTE_NONNUMERIC
import re
import sys

from collections import defaultdict

class defaultdictfromkey(defaultdict):
    def __init__(self, default_factory, unitdict):
        defaultdict.__init__(self, default_factory)
        self._unitdict = unitdict
    def __missing__(self, key):
        unit = self._unitdict.get(key, 'ppt')
        
        var = self.default_factory(key) * factor
        var.units = unit
        return var
        
class flightlogs(PseudoNetCDFFile):
    def __init__(self, pathlike):
        if isinstance(pathlike, str):
            paths = glob(pathlike)
        else:
            paths = pathlike
        paths.sort()
        pfile = self
        pfile._vars = dict()
    
        files = [open(path) for path in paths]
        datas = [np.recfromtxt(f, names = True, case_sensitive = True) for f in files]
        data = np.ma.concatenate(datas)
        desired_unit = dict(O3 = 'ppb', GMAO_TEMP = 'K', PRESS = 'hPa', TEMP = 'K')
        unit_factor = {'ppt': 1e12, 'ppb': 1e9}
        pfile.createDimension('time', data.shape[0])
        for ki, key in enumerate(data.dtype.names):
            typecode = data[key].dtype.char
            if typecode not in ('c', 'S'):
                unit = desired_unit.get(key, 'ppt')
                factor = unit_factor.get(unit, 1)
                values = np.ma.masked_values(data[key], -1000) * factor
            else:
                unit = 'unknown'
                values = data[key]
            pfile.createVariable(key, typecode, dimensions = ('time',), units = unit, values = values)
        
    

if __name__ == '__main__':
    import sys
    bfile1 = flightlogs(sys.argv[1:])
    from matplotlib.mlab import prctile
    for label, key in [('O3', 'O3[:]'), ('NO2', 'NO2[:]')]:
        bvar = eval(key, None, bfile1.variables)
        b2var = eval(key, None, bfile1.variables)
        assert((bvar == b2var).all())
        print >> sys.stdout, '\n%s (BASE: %6.2f)' % (label, bvar.mean())
        print >> sys.stdout, '\n      BASE:',
        prctile(bvar, np.ma.arange(.1, 1., .1)* 100).tofile(sys.stdout, sep = ', ', format = '%6.2f')
    print >> sys.stdout, ''
    