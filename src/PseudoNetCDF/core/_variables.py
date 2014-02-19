from warnings import warn

import numpy as np

class PseudoNetCDFVariable(np.ndarray):
    """
    PseudoNetCDFVariable presents the Scientific.IO.NetCDF.NetCDFVariable interface,
    but unlike that type, provides a contructor for variables that could be used
    without adding it to the parent file
    """
    def __setattr__(self, k, v):
        """
        Set attributes (aka properties) and identify user-defined attributes.
        """
        if not hasattr(self, k) and k[:1] != '_':
            self._ncattrs += (k,)
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

        if hasattr(result, '__dict__'):
            result.__dict__['typecode'] = lambda: typecode
            result.__dict__['dimensions'] = tuple(dimensions)
        else:
            result.__dict__ = {
                'typecode': lambda: typecode,
                'dimensions': tuple(dimensions)
            }

#        object.__setattr__(result, '_ncattrs', ())

        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        if not hasattr(self, '_ncattrs'):
            self._ncattrs = ()
            if obj is None: return
            if hasattr(obj, '_ncattrs'):
                for k in obj._ncattrs:
                    setattr(self, k, getattr(obj, k))
        if not hasattr(self, 'dimensions'):
            if hasattr(obj, 'dimensions'):
                setattr(self, 'dimensions', obj.dimensions)
                self._ncattrs = self._ncattrs[:-1]
        
    
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
    def __setattr__(self, k, v):
        """
        Set attributes (aka properties) and identify user-defined attributes.
        """
        if not hasattr(self, k) and k[:1] != '_':
            self._ncattrs += (k,)
        np.ma.MaskedArray.__setattr__(self, k, v)

    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def ncattrs(self):
        """
        Returns a tuple of attributes that have been user defined
        """
        
        return self._ncattrs

    def swapaxes(self, a1, a2):
        out = np.ma.MaskedArray.swapaxes(self, a1, a2)
        out = self.__array_wrap__(out)
        return out

    def reshape(self, *shape):
        out = np.ma.MaskedArray.reshape(self, *shape)
        out = self.__array_wrap__(out)
        return out
    
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

        result=np.ndarray.view(result, subtype)

        if hasattr(result, '__dict__'):
            result.__dict__['typecode'] = lambda: typecode
            result.__dict__['dimensions'] = tuple(dimensions)
        else:
            result.__dict__ = {
                'typecode': lambda: typecode,
                'dimensions': tuple(dimensions)
            }

#        object.__setattr__(result, '_ncattrs', ())
        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        np.ma.MaskedArray.__array_finalize__(self, obj)
        if not hasattr(self, '_ncattrs'):
            self._ncattrs = ()
            if obj is None: return
            if hasattr(obj, '_ncattrs'):
                for k in obj._ncattrs:
                    setattr(self, k, getattr(obj, k))
        if not hasattr(self, 'dimensions'):
            if hasattr(obj, 'dimensions'):
                setattr(self, 'dimensions', obj.dimensions)
                self._ncattrs = self._ncattrs[:-1]

    def __array_wrap__(self, obj, context = None):
        np.ma.MaskedArray.__array_finalize__(self, obj)
        if not np.isscalar(obj._mask):
            obj._mask = obj._mask.reshape(*obj.shape)
        if hasattr(self, '_ncattrs'):
            for k in self._ncattrs:
                setattr(obj, k, getattr(self, k))
        else:
            self._ncattrs = ()
        
        obj.dimensions = self.dimensions
        return obj
        
    def __getitem__(self, item):
        out = np.ma.MaskedArray.__getitem__(self, item)
        out.dimensions = self.dimensions
        if hasattr(self, '_ncattrs'):
            for k in self._ncattrs:
                setattr(out, k, getattr(self, k))
        else:
            self._ncattrs = ()
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
    
