"""
Written to read files from http://croc.gsfc.nasa.gov/shadoz/
"""
from PseudoNetCDF import PseudoNetCDFFile
import re
import numpy as np
spaces = re.compile(r'\s{2,1000}')
def shadoz(inpath):
    datafile = open(inpath, 'r')
    datalines = datafile.read().split('\n')
    nmeta = int(datalines[0])
    meta = dict([[w.strip() for w in l.split(': ')] for l in datalines[1:nmeta-2]])
    varline, unitline = datalines[nmeta-2:nmeta]
    varnames = spaces.split(varline)
    units = spaces.split(unitline)
    data = np.fromstring('\n'.join(datalines[nmeta:]), sep = ' ').reshape(-1, len(varnames))
    outf = PseudoNetCDFFile()
    outf.createDimension('time', data.shape[0])
    missing = -9999
    for k, v in meta.items():
        setattr(outf, k, v)
        if k == 'Missing or bad values':
            missing = eval(v)
    
    for varname, unit, vals in zip(varnames, units, data.T):
        var = outf.createVariable(varname.replace(' ', '_'), 'f', ('time',))
        var.units = unit
        var.standard_name = varname
        var[:] = np.ma.masked_values(vals[:], missing)
    
    return outf
    
if __name__ == '__main__':
    f = shadoz('/Users/barronh/Downloads/costarica_20160108_V05.1_R.dat')