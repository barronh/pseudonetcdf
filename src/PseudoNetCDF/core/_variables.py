from warnings import warn

import numpy as np

class PseudoNetCDFVariable(np.ndarray):
    """
    PseudoNetCDFVariable presents the Scientific.IO.NetCDF.NetCDFVariable interface,
    but unlike that type, provides a contructor for variables that could be used
    without adding it to the parent file
    """
    __array_priority__ = 100000000.
    def __setattr__(self, k, v):
        """
        Set attributes (aka properties) and identify user-defined attributes.
        """
        if k[:1] != '_' and \
           not k in ('dimensions', 'typecode'):
            if k not in self._ncattrs:
                self._ncattrs += (k, )
        np.ndarray.__setattr__(self, k, v)
    
    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)
    
    def ncattrs(self):
        """
        Returns a tuple of attributes that have been user defined
        """
        
        return self._ncattrs
    
    def __new__(subtype,parent,name,typecode,dimensions,**kwds):
        """
        Creates a variable using the dimensions as defined in
        the parent object

        parent: an object with a dimensions variable
        name: name for variable
        typecode: numpy style typecode
        dimensions: a typle of dimension names to be used from
                    parrent
        kwds: Dictionary of keywords to be added as properties
              to the variable.  **The keyword 'values' is a special
              case that will be used as the starting values of
              the array

        """
        if 'values' in kwds.keys():
            result=kwds.pop('values')
        else:
            shape=[]
            for d in dimensions:
                dim = parent.dimensions[d]

                # Adding support for netCDF3 dimension objects
                if not isinstance(dim, int):
                    dim = len(dim)
                shape.append(dim)

            result=np.zeros(shape,typecode)
        
        result=result[...].view(subtype)

        result.typecode = lambda: typecode
        result.dimensions = tuple(dimensions)
        result._ncattrs = ()
        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        self.typecode = getattr(obj, 'typecode', lambda: self.dtype.char)
        self.dimensions = getattr(obj, 'dimensions', getattr(self, 'dimensions', ()))
        self._ncattrs = getattr(obj, '_ncattrs', getattr(self, '_ncattrs', ()))
        if hasattr(obj, '_ncattrs'):
            for k in obj._ncattrs:
                setattr(self, k, getattr(obj, k))
        
    def swapaxes(self, a1, a2):
        out = np.ndarray.swapaxes(self, a1, a2)
        newdims = list(self.dimensions)
        newdims[a1] = self.dimensions[a2]
        newdims[a2] = self.dimensions[a1]
        out.dimensions = tuple(newdims)
        return out
        
    def getValue(self):
        """
        Return scalar value
        """
        return self.item()

    def assignValue(self,value):
        """
        assign value to scalar variable
        """
        self.itemset(value)

class PseudoNetCDFMaskedVariable(PseudoNetCDFVariable, np.ma.MaskedArray):
    def __new__(subtype,parent,name,typecode = 'f',dimensions = (),**kwds):
        """
        Creates a variable using the dimensions as defined in
        the parent object

        parent: an object with a dimensions variable
        name: name for variable
        typecode: numpy style typecode
        dimensions: a typle of dimension names to be used from
                    parent
        kwds: Dictionary of keywords to be added as properties
              to the variable.  **The keyword 'values' is a special
              case that will be used as the starting values of
              the array

        """
        __array_priority__ = 1000000000.
        if 'values' in kwds.keys():
            result=kwds.pop('values')
        else:
            shape=[]
            for d in dimensions:
                dim = parent.dimensions[d]

                # Adding support for netCDF3 dimension objects
                if not isinstance(dim, int):
                    dim = len(dim)
                shape.append(dim)

            result=np.ma.zeros(shape,typecode)

        result=result.view(subtype)
        result._ncattrs = ()
        result.typecode = lambda: typecode
        result.dimensions = tuple(dimensions)
        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        np.ma.MaskedArray.__array_finalize__(self, obj)
    
    def _update_from(self, obj):
        self.typecode = getattr(obj, 'typecode', lambda: self.dtype.char)
        self.dimensions = getattr(obj, 'dimensions', getattr(self, 'dimensions', ()))
        self._ncattrs = getattr(obj, '_ncattrs', getattr(self, '_ncattrs', ()))
        self._fill_value = getattr(obj, '_fill_value', getattr(self, '_fill_value', -999))
        if hasattr(obj, '_ncattrs'):
            for k in obj._ncattrs:
                setattr(self, k, getattr(obj, k))
        np.ma.MaskedArray._update_from(self, obj)
                
    def __getitem__(self, item):
        out = np.ma.MaskedArray.__getitem__(self, item)
        out._fill_value = self._fill_value
        out = out.view(PseudoNetCDFMaskedVariable)
        if np.isscalar(out): return out
        if hasattr(self, 'dimensions'):
            out.dimensions = self.dimensions
        if hasattr(self, '_ncattrs'):
            for k in self._ncattrs:
                setattr(out, k, getattr(self, k))
        else:
            self._ncattrs = ()
        return out

    def __setattr__(self, k, v):
        """
        Set attributes (aka properties) and identify user-defined attributes.
        """
        if k[:1] != '_' and \
           not k in ('dimensions', 'typecode'):
            if k not in self._ncattrs:
                self._ncattrs += (k, )
        np.ma.MaskedArray.__setattr__(self, k, v)

    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def swapaxes(self, a1, a2):
        out = np.ma.masked_array.swapaxes(self, a1, a2)
        newdims = list(self.dimensions)
        newdims[a1] = self.dimensions[a2]
        newdims[a2] = self.dimensions[a1]
        out.dimensions = tuple(newdims)
        return out
        

    def ncattrs(self):
        """
        Returns a tuple of attributes that have been user defined
        """
        
        return self._ncattrs

    def getValue(self):
        """
        Return scalar value
        """
        return self.item()

    def assignValue(self,value):
        """
        assign value to scalar variable
        """
        self.itemset(value)
  
def PseudoIOAPIVariable(parent,name,typecode,dimensions,**kwds):
    """
    Creates a variable using the dimensions as defined in
    the parent object
    
    parent: an object with a dimensions variable
    name: name for variable
    typecode: numpy style typecode
    dimensions: a typle of dimension names to be used from
                parent
    units: default = none
    long_name: default = name
    var_desc: default = name
    """

    retval = PseudoNetCDFVariable(parent, name, typecode, dimensions, **kwds)

    if not kwds.has_key('units'):
        warn('IOAPI variables must have units; %s has been initialized with "None" units')
        retval.units = 'None'
        
    if not kwds.has_key('long_name'):
        retval.long_name = name.ljust(16)

    if not kwds.has_key('var_desc'):
        retval.var_desc = name.ljust(80)

    return retval
    
