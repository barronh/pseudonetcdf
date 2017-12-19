from ._files import PseudoNetCDFFile

class WrapDict(dict):
    def __init__(self, other):
        self._other = other
    
    def __getitem__(self, k):
        if k in dict.keys(self):
            return dict.__getitem__(self, k)
        elif k in self._other:
            return self._other[k]
        else:
            raise KeyError('{0} not in dictionary or wrapped dictionary'.format(k))
        
    def items(self):
        for k, v in dict.items(self):
            yield k, v
        
        for k, v in self._other.items():
            if not k in dict.keys(self):
                yield k, v
    
    def keys(self):
        for k in dict.keys(self):
            yield k
        
        for k in self._other.keys():
            if not k in dict.keys(self):
                yield k
    def __len__(self):
        return dict.__len__(self) + len(self._other)
    
    def __iter__(self):
        return self.keys()
    
    def __contains__(self, k):
        for myk in self.keys():
            if k == myk: return True
        else:
            return False

class WrapPNC(PseudoNetCDFFile):
    def isMine(path):
        return False
    
    def __init__(self, *args, **kwds):
        #from PseudoNetCDF import pncopen
        self._file = args[0] #pncopen(*args, **kwds)
        for k in self._file.ncattrs():
            setattr(self, k, getattr(self._file, k))
        self.variables = WrapDict(self._file.variables)
        self.dimensions = WrapDict(self._file.dimensions)
    
if __name__ == '__main__':
    from netCDF4 import Dataset
    cpath = 'CONC/CCTM_CONC_v52_hdifupdate_intel17.0_HEMIS_cb6_20160101'
    wfile = WrapPNC(cpath)
