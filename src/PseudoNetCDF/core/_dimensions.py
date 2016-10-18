class PseudoNetCDFDimension(object):
    """
    Dimension object responds like that of netcdf4-python
    """
    def __init__(self, group, name, size):
        self._len = int(size)
        self._unlimited = False
    def isunlimited(self):
        return self._unlimited
    def __len__(self):
        return self._len
    def setunlimited(self, unlimited):
        self._unlimited = unlimited
    
    def __repr__(self):
        out = object.__repr__(self).replace(' at ', ' (len = %d) at' % len(self))
        return out