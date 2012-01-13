from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariable

def jtable(path):
    fobj = file(path)
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
    import re
    header = re.compile('\s+\d+\s+\d+\s+\d+\s*')
    data = []
    headers = []
    for line in lines:
        if header.match(line):
            headers.append(line)
        else:
            data.append(line)
    import numpy as np
    data = np.fromstring(''.join(data), sep = ' ').reshape(nlevels, nlats, nrxns, nangles)
    headers = np.fromstring(''.join(headers), sep = ' ').reshape(nlevels, nlats, nrxns, 3)
    assert((np.diff(headers[..., 0], axis = 0) == 1).all())
    assert((np.diff(headers[..., 1], axis = 1) == 1).all())
    assert((np.diff(headers[..., 2], axis = 2) == 1).all())
    out = PseudoNetCDFFile()
    out.SDATE = jdate
    out.NLAYS = nlevels
    out.NLATS = nlats
    out.NANGLES = nangles
    out.NRXNS = nrxns
    out.createDimension('LAY', nlevels)
    out.variables['LAY'] = PseudoNetCDFVariable(out, 'LAY', 'f', ('LAY',), values = np.array(levels, dtype = 'f'), unit = 'm')
    out.createDimension('LAT', nlats)
    out.variables['LAY'] = PseudoNetCDFVariable(out, 'LAT', 'f', ('LAT',), values = np.array(lats, dtype = 'f'), unit = 'degrees')
    out.createDimension('ANGLE', nangles)
    out.variables['ANGLE'] = PseudoNetCDFVariable(out, 'ANGLE', 'f', ('ANGLE',), values = np.array(angles, dtype = 'f'), unit = 'hours from noon')

    for rxni, rxn in enumerate(rxns):
        out.variables[rxn] = PseudoNetCDFVariable(out, rxn, 'f', ('LAY', 'LAT', 'ANGLE'), values = data[..., rxni, :].copy(), unit = 'hours from noon')
    
    return out

if __name__ == '__main__':
    j = jtable(path = '/Users/barronh/Development/CMAQ/v4.7/data/jproc/Darwin9_i386/JTABLE_2001203')