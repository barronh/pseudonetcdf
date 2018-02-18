from warnings import warn
from collections import OrderedDict

import numpy as np

class PseudoNetCDFVariable(np.ndarray):
    """
    PseudoNetCDFVariable presents the Scientific.IO.NetCDF.NetCDFVariable interface,
    but unlike that type, provides a contructor for variables that could be used
    without adding it to the parent file
    """
    __array_priority__ = 10000000.
    def xarray(self, iscoord = False):
        """
Experimental function
Returns:
    out - xarray.DataArray object with dimensions and coordinates
        dims - are set by self.get_coords()
        coords - are set by self.get_coords()

Notes:
    When data is 2-d, the dimensions
        """
        import xarray as xr
        myname = self._name
        if iscoord:
            coordd = {}
        else:
            coordd = dict([(k, self.get_coord(k)) for k in self.get_coord_names() if k != myname])
        attrs = dict([(k, getattr(self, k)) for k in self.ncattrs()])
        return xr.DataArray(self[:], dims = self.dimensions, coords = coordd, attrs = attrs)

    def get_coord_names(self):
        dims = list(self.dimensions)
        ndims = len(dims)
        coords = getattr(self, 'coordinates', '').split()
        ncoords = len(coords)
        maxdims = max(ndims, ncoords);
        coords = ((ndims * [None]) + coords)[-maxdims:]
        dims = (dims + ncoords * [None])[:maxdims]
        coords = [dimn if coordn is None else coordn for dimn, coordn in zip(dims, coords)]
        return coords
    
    def get_coord(self, coordn):
        if coordn in self._parent.variables:
            v = self._parent.variables[coordn]
            if len(v.dimensions) > 1 and coordn in v.dimensions:
                warn(coordn + ' has dimensions ' + str(v.dimensions) + '; so it will use a standard linear coordinate. Consider manually adding coords')
                coordi = list(self.dimensions).index(coordn)
                coordv = np.arange(self.shape[coordi])
            else:
                coordv = (v.dimensions, v.xarray(iscoord = True))
        else:
            coordi = list(self.dimensions).index(coordn)
            coordv = np.arange(self.shape[coordi])
        return coordv
     
    def get_coords(self):
        coords = self.get_coord_names()
        coordsd = {}
        for coordn in coords:
             coordsd[coordn] = self.get_coord(coordn)
        return coordsd

    def __repr__(self):
        out = str(self)
        return object.__repr__(self).replace(' at ', '\n' + out + ' at ')
    
    def __str__(self):
        namekeys = ['_name', 'name', 'standard_name', 'long_name']
        for nck in self.ncattrs():
            if 'name' in nck:
                namekeys.append(nck)
        
        for nck in namekeys:
            if hasattr(self, nck):
                var_name = getattr(self, nck)
                break
        else:
            var_name = 'unknown'
                                
        var_type = dict(float32='float', \
                        float64='double', \
                        int32='integer', \
                        uint32='integer', \
                        int64='long', \
                        bool='bool', \
                        string8='char', \
                        string80='char').get(self.dtype.name, self.dtype.name)
        out = ""
        indent = "    "
        startindent = ""
        out += startindent + 1*indent+("%s %s%s; // shape: %s\n" % (var_type, var_name,str(self.dimensions).replace('u\'', '').replace('\'','').replace(',)',')'), str(self.shape)))
        for prop_name in self.ncattrs():
            prop = getattr(self, prop_name)
            out += startindent + 2*indent+("%s:%s = %s ;\n" % (var_name,prop_name,repr(prop).replace("'", '"')))
        
        out += 'array: '
        out += self.array().__str__()
        return out

    def array(self):
        """
        Return parent type view object
        """
        return self.view(type = np.ndarray)
    
    def __setattr__(self, k, v):
        """
        Set attributes (aka properties) and identify user-defined attributes.
        """
        if k[:1] != '_' and \
           not k in ('dimensions', 'typecode'):
            if k not in self._ncattrs:
                self._ncattrs += (k, )
        object.__setattr__(self, k, v)
    
    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def setncattr(self, k, v):
        return setattr(self, k, v)
    
    def getncattr(self, k):
        return getattr(self, k)
    
    def setncatts(self, attdict):
        for pk, pv in attdict.items():
            setattr(self, pk, pv)

    def getncatts(self):
        outd = OrderedDict()
        for pk in self.ncattrs():
            outd[pk] = getattr(self, pk)
        return outd
    
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
                    parent
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
        result._parent= parent
        result._name = name
        for k,v in kwds.items():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        if obj is None: return
        _parent = getattr(obj, '_parent', object())
        object.__setattr__(self, '_parent', _parent)
        _name = getattr(obj, '_name', 'unknown')
        object.__setattr__(self, '_name', _name)
        ntypecode = getattr(obj, 'typecode', lambda: self.dtype.char)
        object.__setattr__(self, 'typecode', ntypecode)
        ndimensions = getattr(obj, 'dimensions', lambda: self.dtype.char)
        object.__setattr__(self, 'dimensions', ndimensions)
        nncattrs = getattr(obj, '_ncattrs', getattr(self, '_ncattrs', ()))
        object.__setattr__(self, '_ncattrs', nncattrs) 
        if hasattr(obj, '_ncattrs'):
            for k in nncattrs:
                if not hasattr(self, k):
                    object.__setattr__(self, k, getattr(obj, k))
        
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

            result=np.ma.zeros(shape, dtype = 'S1' if typecode == 'c' else typecode, fill_value = kwds.get('fill_value', None))

        result=result.view(subtype)
        result._ncattrs = ()
        result._parent = parent
        result._name = name
        result.typecode = lambda: typecode
        result.dimensions = tuple(dimensions)
        for k,v in kwds.items():
            setattr(result,k,v)
        return result

    def array(self):
        """
        Return parent type view object
        """
        return self.view(type = np.ma.masked_array)
    
    def __array_finalize__(self, obj):
        np.ma.MaskedArray.__array_finalize__(self, obj)
    
    def _update_from(self, obj):
        dt = self.dtype.char
        self.typecode = getattr(obj, 'typecode', lambda: ('c' if dt == 'S' else dt))
        self.dimensions = getattr(obj, 'dimensions', getattr(self, 'dimensions', ()))
        self._ncattrs = getattr(obj, '_ncattrs', getattr(self, '_ncattrs', ()))
        self._fill_value = getattr(obj, '_fill_value', getattr(self, '_fill_value', -999))
        if hasattr(obj, '_ncattrs'):
            for k in obj._ncattrs:
                if k in ('fill_value',): continue
                setattr(self, k, getattr(obj, k))
        np.ma.MaskedArray._update_from(self, obj)
                
    def __getitem__(self, item):
        out = np.ma.MaskedArray.__getitem__(self, item)
        try: out._fill_value = self._fill_value
        except: pass
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

    if not 'units' in kwds:
        warn('IOAPI variables must have units; %s has been initialized with "None" units')
        retval.units = 'None'
        
    if not 'long_name' in kwds:
        retval.long_name = name.ljust(16)

    if not 'var_desc' in kwds:
        retval.var_desc = name.ljust(80)

    return retval
    
