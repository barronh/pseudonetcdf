from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable
import re
import numpy as np

class jtable(PseudoNetCDFFile):
    def __init__(self, path):
        self.dimensions = {}
        self.variables = {}
        fobj = open(path)
        jdate = int(fobj.readline().split()[0])
        nlevels = int(fobj.readline().split()[0])
        levels = [float(x) for x in fobj.readline().split()]
        nlats = int(fobj.readline().split()[0])
        lats = [float(x) for x in fobj.readline().split()]
        nangles = int(fobj.readline().split()[0])
        angles = [float(x) for x in fobj.readline().split()]
        nrxns = int(fobj.readline().split()[0])
        scaling = dict()
        rxns = []
        for ri in range(nrxns):
            tname, scale = fobj.readline().split(',')
            tname = tname.replace("'", "").strip()
            scaling[tname] = eval(scale.strip())
            rxns.append(tname)
        lines = fobj.readlines()
        header = re.compile('\s+\d+\s+\d+\s+\d+\s*')
        data = []
        headers = []
        for line in lines:
            if header.match(line):
                headers.append(line)
            else:
                data.append(line)
        data = np.fromstring(''.join(data), sep = ' ').reshape(nlevels, nlats, nrxns, nangles)
        headers = np.fromstring(''.join(headers), sep = ' ').reshape(nlevels, nlats, nrxns, 3)
        assert((np.diff(headers[..., 0], axis = 0) == 1).all())
        assert((np.diff(headers[..., 1], axis = 1) == 1).all())
        assert((np.diff(headers[..., 2], axis = 2) == 1).all())
        self.SDATE = jdate
        self.NLAYS = nlevels
        self.NLATS = nlats
        self.NANGLES = nangles
        self.NRXNS = nrxns
        self.createDimension('LAY', nlevels)
        self.variables['LAY'] = PseudoNetCDFVariable(self, 'LAY', 'f', ('LAY',), values = np.array(levels, dtype = 'f'), units = 'm')
        self.createDimension('LAT', nlats)
        self.variables['LAT'] = PseudoNetCDFVariable(self, 'LAT', 'f', ('LAT',), values = np.array(lats, dtype = 'f'), units = 'degrees')
        self.createDimension('ANGLE', nangles)
        self.variables['ANGLE'] = PseudoNetCDFVariable(self, 'ANGLE', 'f', ('ANGLE',), values = np.array(angles, dtype = 'f'), units = 'hours from noon')
    
        for rxni, rxn in enumerate(rxns):
            self.variables[rxn] = PseudoNetCDFVariable(self, rxn, 'f', ('LAY', 'LAT', 'ANGLE'), values = data[..., rxni, :].copy(), units = 's**-1')

    def interpj(self, key_or_var, lay, lat, angle):
        from scipy.interpolate import griddata
        lays = self.variables['LAY']
        lats = self.variables['LAT']
        angles = self.variables['ANGLE']
        points = np.array([[[(lay_, lat_, angle_) for angle_ in angles] for lat_ in lats] for lay_ in lays]).reshape(-1, 3)
        if isinstance(key_or_var, str):
            v = self.variables[key_or_var].ravel()
        else:
            v = np.array(key_or_var)
            assert(key_or_var.shape == (self.NLAYS, self.NLATS, self.NANGLES))
        if isinstance(lay, (int, float)):
            idx = np.array((lay, lat, angle), ndmin = 2)
        else:
            idx = np.array((lay, lat, angle)).swapaxes(0,1)
            assert(idx.shape[0] == 2)
        return griddata(points, v, idx)[0]
        

if __name__ == '__main__':
    j = jtable(path = '/Volumes/LaCie/JTABLE_1985172')
    v = j.variables['HCHO_M_SAPRC99']
    np.set_printoptions(precision = 5)
    print v[:2,:2,:2]
    print j.interpj('HCHO_M_SAPRC99', 0,10,0)
    print j.interpj('HCHO_M_SAPRC99', 0,10,1)
    print j.interpj('HCHO_M_SAPRC99', 0,20,1)
    print j.interpj('HCHO_M_SAPRC99', 1000,20,1)
    print j.interpj('HCHO_M_SAPRC99', 0,10,.5)
    print j.interpj('HCHO_M_SAPRC99', 0,15,.5)
    print j.interpj('HCHO_M_SAPRC99', 500,15,.5)