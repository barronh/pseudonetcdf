__all__ = ['bpch']

import os

# part of the default Python distribution
from collections import defaultdict
from warnings import warn

# numpy is a very common installed library
from numpy import fromfile, memmap, dtype, arange, zeros

# PseudoNetCDF is my own home grown
# https://dawes.sph.unc.edu/trac/PseudoNetCDF
from PseudoNetCDF import PseudoNetCDFVariable, PseudoNetCDFFile


# These variables define the binary format of the header blocks
# and are only for internal
_general_header_type = dtype('>i4, S40, >i4, >i4, S80, >i4')
_datablock_header_type = dtype('>i4, S20, 2>f4, >i4, >i4, >i4, >i4, S40, >i4, S40, >f8, >f8, S40, 3>i4, 3>i4, >i4, >i4')
_first_header_size = _general_header_type.itemsize + _datablock_header_type.itemsize

class defaultdictfromkey(defaultdict):
    """
    defaultdictfromkey dynamically produces dictionary items
    for keys, using the default_factor function called
    with the key as an argument
    """

    def __missing__(self, key):
        """
        __missing__(key) # Called by __getitem__ for missing key; pseudo-code:
        if self.default_factory is None: raise KeyError((key,))
        self[key] = value = self.default_factory()
        return value
        """
        return self.default_factory(key)

class defaultdictfromthesekeys(defaultdict):
    """
    defaultdictfromthesekeys dynamically produces dictionary items
    for known keys, using the default_factor function called
    with the key as an argument
    """
    def __init__(self, keys, default_factory = None):
        """
        keys - iterable of keys that default_factory can produce
        default_factory - function that takes a key as an argument
                          to create a dictionary item
        """
        self._keys = set([k for k in keys])
        defaultdict.__init__(self, default_factory)
    
    def __iter__(self):
        for i in self._keys:
            yield i
            
    def iterkeys(self):
        for i in self._keys:
            yield i
    
    def itervalues(self):
        for k in self.iterkeys():
            yield self[k]
    
    def iteritems(self):
        for k in self.iterkeys():
            yield (k, self[k])
    
    def keys(self):
        return [k for k in self]
    
    def __setitem__(self, key, value):
        val = defaultdict.__setitem__(self, key, value)
        self._keys.add(key)
        return val
    
    def __delitem__(self, key):
        self._keys.discard(key)
        return defaultdict.__del__(self, key)
    
    def pop(self, key):
        val = defaultdict.pop(self, key)
        self._keys.discard(key)
        return val
        
    def __missing__(self, key):
        """
        __missing__(key) # Called by __getitem__ for missing key; pseudo-code:
        if self.default_factory is None: raise KeyError((key,))
        self[key] = value = self.default_factory()
        return value
        """
        if key in self._keys:
            return self.default_factory(key)
        else:
            raise KeyError("%s not found in %s" % (key, ))

    for k in '__setitem__ __delitem__ pop __iter__ iterkeys itervalues iteritems keys'.split():
        exec('%s.__doc__ = defaultdict.%s.__doc__' % (k, k))

class defaultpseudonetcdfvariable(defaultdictfromthesekeys):
    """
    Overwrites __repr__ function to show variables
    """
    def __repr__(self):
        out = "{"
        for k in self._keys:
            out += "\n  '%s': PseudoNetCDFVariable(...)" % k
        out += "\n}"
        return out

    def __str__(self):
        return self.__repr__()

class _diag_group(PseudoNetCDFFile):
    """
    This object acts as a PseudoNetCDF file that gets data from the parent object.
    """
    def __init__(self, parent, groupname, groupvariables):
        """
        parent - a PseudoNetCDFFile
        groupname - a string describing the group
        """
        template = '%s_%%s' % groupname
        def getvar(key):
            try:
                return parent.variables[template % key]
            except (KeyError, ValueError):
                return parent.variables[key]
        self._parent = parent
        self.variables = defaultpseudonetcdfvariable(list(groupvariables) + ['tau0', 'tau1', 'LAT', 'LON'], getvar)
    
    def __getattr__(self, key):
        try:
            return object.__getattr__(self, key)
        except AttributeError:
            return getattr(self._parent, key)
            
# This class is designed to operate like a dictionary, but
# dynamically create variables to return to the user
class _tracer_lookup(defaultpseudonetcdfvariable):
    """
    _tracer_lookup: finds geos_chem tracer indices from names and returns
                    netcdf like variable
    """
    def __init__(self, parent, datamap, tracerinfo, diaginfo, keys):
        """
        parent: NetCDF-like object to serve dimensions
        datamap: array of pre-dimensioned and datatyped values dim(tstep, i, j, k)
        tracerinfo: dictionary of tracer data keyed by ordinal
        keys: list of keys to serve
        """
        self._tracer_data = tracerinfo
        self._diag_data = diaginfo
        self._memmap = datamap
        self._parent = parent
        self._keys = keys + ['tau0', 'tau1', 'LAT', 'LON']
        
    def __missing__(self, key):
        if key == 'LAT':
            j = arange(len(self._parent.dimensions['J']) + 1)
            data = j * self._parent.modelres[1] - 90
            if self._parent.halfpolar == 1:
                data[1:] -= 1
                data[-1] -= 1
            dims = ('J',)
            dtype = 'i'
            units = 'degrees north'
        elif key == 'LON':
            i = arange(len(self._parent.dimensions['I']) + 1)
            xres = self._parent.modelres[0]
            data = i * xres - (180 + xres / 2.) * self._parent.center180
            dims = ('I',)
            dtype = 'i'
            units = 'degrees east'
        elif key == 'tau0':
            tmp_key = self._keys[0]
            data = self._memmap[tmp_key]['header']['f10']
            dims = ('T',)
            dtype = 'i'
            units = 'hours since 0 GMT 1/1/1985'
        elif key == 'tau1':
            tmp_key = self._keys[0]
            data = self._memmap[tmp_key]['header']['f11']
            dims = ('T',)
            dtype = 'i'
            units = 'hours since 0 GMT 1/1/1985'
        else:
            dims = ('T', 'K', 'J', 'I')
            dtype = 'f'
            header = self._memmap[key]['header'][0]
            sl, sj, si = header['f14'][::-1] - 1
            group = header['f1']
            offset = self._diag_data[group]['offset']
            ord = header['f8'] + offset
            base_units = header['f9']
            scale = self._tracer_data[ord]['SCALE']
            carbon = self._tracer_data[ord]['C']
            units = self._tracer_data[ord]['UNIT']
            if carbon != 1 and 'C' not in units:
                warn("Scaling %s by carbon, but unit does not indicate carbon" % key)
            tmp_data = self._memmap[key]['data']
            assert((tmp_data['f0'] == tmp_data['f2']).all())
            data = tmp_data['f1'] * scale * carbon
            if any([sl != 0, sj != 0, si != 0]):
                warn("%s is a subset variable" % key)
                nl, nj, ni = header['f14'][::-1]
                tmp_data = zeros((data.shape[0], nl, nj, ni), dtype = data.dtype)
                el, ej, ei = data.shape[1:]
                el += sl
                ej += sj
                ei += si
                tmp_data[:, sl:el, sj:ej, si:ei] = data[:]
                data = tmp_data
        return PseudoNetCDFVariable(self._parent, key, dtype, dims, units = units, long_name = key, var_desc = key, values = data)
            
    def __del__(self):
        del self._memmap


class bpch(PseudoNetCDFFile):
    """
    NetCDF-like class to interface with GEOS-Chem binary punch files
    
    f = bpch(path_to_binary_file)
    dim = f.dimensions[dkey] # e.g., dkey = 'I'
    
    # There are two ways to get variables.  Directly from
    # the file using the long name
    var = f.variables[vkey] # e.g., vkey = 'IJ-AVG-$_NOx'
    
    # Or through a group using the short name
    g = f.groups[gkey] # e.g., gkey = 'IJ-AVG-$'
    var = g.variables[vkey] # e.g., vkey = 'NOx'

    # The variable returned is the same either way    
    print f.dimensions
    print var.unit
    print var.dimensions
    print var.shape
    
    """
    def __init__(self, bpch_path, tracerinfo = None, diaginfo = None, mode = 'r'):
        """
        bpch_path: path to binary punch file
        tracerinfo: path to ascii file with tracer definitions
        diaginfo: path to ascii file with diagnostic group definitions
        mode : {'r+', 'r', 'w+', 'c'}, optional
         |      The file is opened in this mode:
         |  
         |      +------+-------------------------------------------------------------+
         |      | 'r'  | Open existing file for reading only.                        |
         |      +------+-------------------------------------------------------------+
         |      | 'r+' | Open existing file for reading and writing.                 |
         |      +------+-------------------------------------------------------------+
         |      | 'w+' | Create or overwrite existing file for reading and writing.  |
         |      +------+-------------------------------------------------------------+
         |      | 'c'  | Copy-on-write: assignments affect data in memory, but       |
         |      |      | changes are not saved to disk.  The file on disk is         |
         |      |      | read-only.                                                  |
         |      +------+-------------------------------------------------------------+        
        """
        # Read binary data for general header and first datablock header
        header_block = fromfile(bpch_path, dtype = 'bool', count = _first_header_size)
        
        # combine data for convenience
        header = tuple(header_block[:_general_header_type.itemsize].view(_general_header_type)[0]) + \
                 tuple(header_block[_general_header_type.itemsize:].view(_datablock_header_type)[0])
        
        # Verify that all Fortran unformatted buffers match 
        try:
            assert(header[0] == header[2])
            assert(header[3] == header[5])
        except AssertionError:
            raise ValueError("BPCH Files fails header check")
        
        # Assign data from header to global attributes
        self.ftype = header[1]
        self.toptitle = header[4]
        self.modelname, self.modelres, self.halfpolar, self.center180 = header[7:11]
        dummy, dummy, dummy, self.start_tau0, self.start_tau1, dummy, dim, dummy, dummy = header[13:-1]
        self.dimensions = dict(zip('I J K'.split(), dim))
        
        tracerinfo = tracerinfo or os.path.join(os.path.dirname(bpch_path), 'tracerinfo.dat')
        if os.path.exists(tracerinfo):
            tracer_data = dict([(int(l[52:61].strip()), dict(NAME = l[:8].strip(), FULLNAME = l[9:39].strip(), MOLWT = float(l[39:49]), C = int(l[49:52]), TRACER = int(l[52:61]), SCALE = float(l[61:71]), UNIT = l[72:].strip())) for l in file(tracerinfo).readlines() if l[0] not in ('#', ' ')])
            tracer_names = dict([(k, v['NAME']) for k, v in tracer_data.iteritems()])
            add_unit = False
        else:
            warn('Reading file without tracerinfo.dat means that names and scaling are unknown')
            tracer_data = defaultdict(lambda: dict(SCALE = 1., C = 1.))
            tracer_names = defaultdictfromkey(lambda key: key)
            add_unit = True
        
        diaginfo = diaginfo or os.path.join(os.path.dirname(bpch_path), 'diaginfo.dat')
        if os.path.exists(diaginfo):
            diag_data = dict([(l[9:49].strip(), dict(offset = int(l[:8]), desc = l[50:].strip())) for l in file(diaginfo).read().strip().split('\n') if l[0] != '#'])
        else:
            warn('Reading file without diaginfo.dat loses descriptive information')
            diag_data = defaultdictfromkey(lambda key: dict(offset = 0, desc = key))
            
        if len(tracer_names) == 0 and not add_unit:
            raise IOError("Error parsing %s for Tracer data")
        file_size = os.stat(bpch_path).st_size
        offset = _general_header_type.itemsize
        data_types = []
        first_header = None
        keys = []
        self._groups = defaultdict(set)
        while first_header is None or \
              offset < file_size:
            header = memmap(bpch_path, offset = offset, shape = (1,), dtype = _datablock_header_type, mode = mode)[0]
            
            group = header[7].strip()
            goffset = diag_data[group]['offset']
            tracer_number = header[8]
            unit = header[9].strip()
            if add_unit:
                tracer_data[tracer_number]['UNIT'] = unit
            
            try:
                tracername = tracer_names[tracer_number + goffset]
            except:
                # There are some cases like with the adjoint where the tracer does not have
                # a tracerinfo.dat line.  In this case, the name matches the tracer with that
                # number (no offset).  The scaling, however, is not intended to be used.
                # The unit for adjoint, for instance, is unitless.
                tracername = tracer_names[tracer_number]
                tracer_data[tracer_number] = dict(SCALE = 1., C = tracer_data[tracer_number + goffset]['C'], UNIT = unit)

            self._groups[group].add(tracername)
            if first_header is None:
                first_header = header
            elif (header[7], header[8]) == (first_header[7], first_header[8]):
                break
            dim = header[13][::-1]
            start = header[14][::-1]
            data_type = dtype('>i4, %s>f4, >i4' % str(tuple(dim[:])))
            assert(data_type.itemsize == header[-2])
            data_types.append(data_type)
            keys.append('%s_%s' % (group, tracername))
            offset += _datablock_header_type.itemsize + header[-2]

        time_type = dtype([(k, dtype([('header', _datablock_header_type), ('data', d)])) for k, d in zip(keys, data_types)])
        # load all data blocks
        datamap = memmap(bpch_path, dtype = time_type, offset = _general_header_type.itemsize, mode = mode)
        
        # Create variables and dimensions
        self.variables = _tracer_lookup(parent = self, datamap = datamap, tracerinfo = tracer_data, diaginfo = diag_data, keys = keys)
        self.createDimension('T', self.variables['tau0'].shape[0])
        self.groups = dict([(k, _diag_group(self, k, v)) for k, v in self._groups.iteritems()])


    def __del__(self):
        del self.variables
      
    def __repr__(self):
        return PseudoNetCDFFile.__repr__(self) + str(self.variables)

if __name__ == '__main__':
    from numpy import median, indices, arange, meshgrid

    # Example: file open and variable aquisition
    path_to_test_file = 'restart.geos5.2005010100'
    path_to_test_file = 'ctm.bpch'
    f = None
    while f is None:
        try:
            print path_to_test_file
            f = bpch(path_to_test_file)
        except:
            path_to_test_file = raw_input('Enter path to a valid GEOS-Chem file')
    
    group_key = f.groups.keys()[0]
    g = f.groups[group_key]
    i = 0; var_key = 'LAT'
    while var_key in 'LAT LON tau0 tau1'.split():
        var_key = g.variables.keys()[i]
        i += 1
    var = g.variables[var_key]
    
    # Example: variable metadata print
    print var.long_name
    print var.dimensions
    print var.shape
    
    # Example: time/layer slice
    layer1 = var[0, 0, :, :]
    
    # Example: data statistics
    print layer1.min(), layer1.mean(), median(layer1), layer1.max(), var.units.strip()
    try:
        # Example: spatial plotting
        from matplotlib import use
        use('Agg')
        from pylab import figure, xticks, yticks, title, colorbar
        from mpl_toolkits.basemap import Basemap
        
        # I'm actually not sure what the right projection for this 
        # data is.  So the next few lines might need to change.
        m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-180 - f.modelres[1] * f.center180 / 2.,\
                    urcrnrlon=180 - f.modelres[1] * f.center180 / 2.,\
                    resolution='c')
    
        lat = f.variables['LAT']
        lon = f.variables['LON']
        x, y = meshgrid(*m(lon, lat))
        
        fig = figure(figsize = (9,4))
        ax = fig.add_axes([.05, .1, .9, .8])
        m.drawcountries()
        m.drawstates()
        m.drawcoastlines()
        parallels = arange(-90,91,15)
        meridians = arange(-180,180,30)
        m.drawparallels(parallels)
        m.drawmeridians(meridians)
        poly = m.pcolor(lon, lat, layer1)
        cb = colorbar(poly, ax = m.ax)
        cb.ax.set_xlabel(var.units.strip())
        xticks(meridians)
        yticks(parallels)
        title('%s (Projection might be wrong)' % var_key)
        fig_path = 'layer1_%s.png' % var_key
        fig.savefig(fig_path)
        print("Examine test figure %s" % fig_path)
    except Exception, e:
        print("Unable to produce test figure (maybe you don't have matplotlib or basemap); " + str(e))