from ._files import PseudoNetCDFFile
from collections import OrderedDict

class WrapDict(OrderedDict):
    def __init__(self, other = {}):
        self._other = other
        self._mine = OrderedDict()
    
    def __getitem__(self, k):
        if k in self._mine:
            return self._mine[k]
        elif k in self._other:
            return self._other[k]
        else:
            raise KeyError('{0} not in dictionary or wrapped dictionary'.format(k))
        
    def items(self):
        for k, v in self._other.items():
            if not k in self._mine:
                yield k, v
        
        for k, v in self._mine.items():
            yield k, v
        
    def keys(self):
        for k in self._other.keys():
            if not k in self._mine:
                yield k
        
        for k in self._mine.keys():
            yield k
        
    def __len__(self):
        return len(self._mine) + len(self._other)
    
    def __iter__(self):
        return self.keys()
    
    def __contains__(self, k):
        for myk in self.keys():
            if k == myk: return True
        else:
            return False
    
    def __setitem__(self, k, v):
        self._mine[k] = v

    def __delitem__(self, k):
        if k in self._mine:
            del self._mine[k]
        
        if k in self._other:
            del self._other[k]
    
    def __repr__(self):
        outf = OrderedDict()
        for key in self.keys():
            outf[key] = self[key]
        return outf.__repr__()
    
    def __str__(self):
        outf = OrderedDict()
        for key in self.keys():
            outf[key] = self[key]
        return outf.__str__()
    
class WrapPNC(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, path):
        return False
    
    def __init__(self, *args, **kwds):
        #from PseudoNetCDF import pncopen
        self._file = args[0] #pncopen(*args, **kwds)
        for k in self._file.ncattrs():
            setattr(self, k, getattr(self._file, k))
        self.dimensions = WrapDict(self._file.dimensions)
        self.variables = WrapDict(self._file.variables)
    
if __name__ == '__main__':
    from netCDF4 import Dataset
    cpath = 'CONC/CCTM_CONC_v52_hdifupdate_intel17.0_HEMIS_cb6_20160101'
    wfile = WrapPNC(cpath)
