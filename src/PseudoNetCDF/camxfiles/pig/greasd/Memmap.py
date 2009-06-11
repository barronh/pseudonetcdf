from numpy import memmap, dtype
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
import os

class greasd(PseudoNetCDFFile):
    def __init__(self, pig_path, mode = 'r'):
        file_size = os.path.getsize(pig_path)
        tmpmap = memmap(pig_path, dtype = '>i', shape = 7)
        self.SDATE = tmpmap[1]
        self.STIME = tmpmap[2:3].view('>f')[0]
        self.NPIG = tmpmap[3]
        record_length = int(tmpmap[-1]/4)
        per_pig = record_length / self.NPIG
        
        pig_dtype = dtype(
                         dict(
                           names = 'SPAD1 IDATE TIME NPIG IDUM EPAD1'.split()+['SPAD2']+['P%d_%d' % (i / per_pig,i % per_pig) for i in range(record_length)]+['EPAD2'], 
                           formats = '>i >i >f >i >i >i >i'.split() + ['>f'] * record_length + ['>i']
                         )
                       )
        items = float(file_size) / float(pig_dtype.itemsize)
        assert(items == int(items))
        self.NSTEPS = items = int(items)
        self.__memmap__ = memmap(pig_path, shape = items, dtype = pig_dtype)
        self.dimensions = dict(TSTEP = items, NPIG = self.NPIG)
        self.variables = PseudoNetCDFVariables(self.__variables, pig_dtype.fields.keys())
    def __variables(self,k):
        return PseudoNetCDFVariable(self, k, 'f', ('TSTEP',), values = self.__memmap__[k], var_desc = k.ljust(16), long_name = k.ljust(16))