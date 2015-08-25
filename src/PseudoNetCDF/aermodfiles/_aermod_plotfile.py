import numpy as np
from PseudoNetCDF.sci_var import PseudoNetCDFFile
names = 'X,Y,CONC,ZELEV,ZHILL,ZFLAG,dum1,AVE,dum2,GRP,DATE,dum3,NETID'.split(',')
units = 'm,m,micrograms/m**3,m,m,m,N/A,N/A,YYMMDDHH,N/A'.split(',')
#(3(1X,F13.5),3(1X,F8.2),2X,A6,2X,A8,2X,I8.8,2X,A8) 
delimiter = [14] * 3 + [9] * 3 + [2, 6, 2, 8, 10, 2, 8]
StrLen = 10
class reader(PseudoNetCDFFile):
    def isMine(self, path):
        return open(path).read(len('* AERMOD')) == '* AERMOD'
    def __init__(self, path):
        #import pdb; pdb.set_trace()
        self._data = np.recfromtxt(path, names = names, delimiter = delimiter, comments = '*')
        self.createDimension('receptor', self._data.size)
        self.createDimension('StrLen', StrLen)
        for k, u in zip(names, units):
            if k in 'dum1 dum2 dum3'.split(): continue
            vals = self._data[k]
            dt = vals.dtype.char
            dims = ('receptor',)
            if dt == 'S':
                vals = np.char.ljust(np.char.strip(vals), StrLen).view('S1').reshape(-1, StrLen)
                dims = ('receptor', 'StrLen')
                dt = 'S1'
            var = self.createVariable(k, dt, dims)
            var.units = u
            var[:] = vals
    
        