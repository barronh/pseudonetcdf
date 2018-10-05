__all__ = ['PseudoNetCDFDimension']


class PseudoNetCDFDimension(object):
    """
    Dimension object responds like that of netcdf4-python
    """

    def __init__(self, group, name, size):
        self._len = int(size)
        self._unlimited = False
        self._name = name

    def isunlimited(self):
        """
        Returns
        -------
        out : boolean
            is the dimension unlimited
        """
        return self._unlimited

    def __len__(self):
        return self._len

    def setunlimited(self, unlimited):
        """
        Parameters
        ----------
        unlimited : boolean
            largely decorative in PseudoNetCDF, but important in persisting
        """
        self._unlimited = unlimited

    def __repr__(self):
        out = object.__repr__(self).replace(
            ' at ', ' (name = %s, len = %d) at' % (self._name, len(self)))
        return out
