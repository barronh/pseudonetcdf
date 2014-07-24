import numpy as np
import sys
from PseudoNetCDF import PseudoNetCDFFile
def raob(inpath):
    nfields, = np.fromfile(inpath, dtype = '>i4', count = 1) / 4.
    fblock = np.fromfile(inpath, dtype = '>i4')
    block = fblock[1:-2]
    block = block.reshape(-1, nfields + 1)[:, :-1];
    names, units = np.char.strip(block[:2].view('S4'));
    datas = (np.ma.masked_values(block[3, :], -2139062144).filled(-999) / 10.**block[2]);
    data = (np.ma.masked_values(block[4:, :], -2139062144).filled(-999) / 10.**block[2]);
    outf = PseudoNetCDFFile()
    outf.createDimension('level', data.shape[0])
    for k, u, v, vs in zip(names, units, data.T, datas):
        outf.createVariable(k, 'f', ('level',), units = u, values = v)
        outf.createVariable(k + '_FIRST', 'f', ('level',), units = u, values = vs)
        
    #np.savetxt(sys.stdout, (names, units), delimiter = ', ', fmt = '%7s')
    #np.savetxt(sys.stdout, data, delimiter = ', ', fmt = '%7.1f')
    return outf

if __name__ == '__main__':
    ifile0 = raob(sys.argv[1])