__all__ = ['bpch']

# collections is part of the default Python distribution
from collections import defaultdict
import os
from warnings import warn

# numpy is a very common installed library
from numpy import fromfile, memmap, dtype, median, arange

# PseudoNetCDF is my own home grown
# https://dawes.sph.unc.edu/trac/PseudoNetCDF
from PseudoNetCDF import PseudoNetCDFVariable, PseudoNetCDFFile

# These variables define the binary format of the header blocks
# and are only for internal
_general_header_type = dtype('>i4, S40, >i4, >i4, S80, >i4')
_datablock_header_type = dtype('>i4, S20, 2>f4, >i4, >i4, >i4, >i4, S40, >i4, S40, >f8, >f8, S40, 3>i4, 3>i4, >i4, >i4')
_first_header_size = _general_header_type.itemsize + _datablock_header_type.itemsize

# This class is designed to operate like a dictionary, but
# dynamically create variables to return to the user
class _tracer_lookup(defaultdict):
    """
    _tracer_lookup: finds geos_chem tracer indices from names and returns
                    netcdf like variable
    """
    def __init__(self, parent, headers, data, tracerinfo):
        """
        parent: NetCDF-like object to serve dimensions
        header: array of datablock headers
        data: array of pre-dimensioned and datatyped values dim(var, i, j, k)
        tracerinfo: path to file with tracer definitions
        """
        try:
            self._tracer_data = dict([(int(l[52:61].strip()), dict(NAME = l[:8].strip(), FULLNAME = l[9:39].strip(), MOLWT = float(l[39:49]), C = int(l[49:52]), TRACER = int(l[52:61]), SCALE = float(l[61:71]), UNIT = l[72:].strip())) for l in file(tracerinfo).readlines() if l[0] not in ('#', ' ')])
            self._tracer_names = dict([(k, v['NAME']) for k, v in self._tracer_data.iteritems()])
        except:
            self._tracer_names = {}
        
        if len(self._tracer_names) == 0:
            raise IOError, "Error parsing %s for Tracer data"
            
        self._headers = headers
        self._data = data
        self._parent = parent
        self._prune()
        self._keys += ['tau0', 'tau1', 'LAT', 'LON']

    def keys(self):
        return [i for i in self._keys]
        
    def iteritems(self):
        return [(i, self[i]) for i in self._keys]

    def values(self):
        return [self[i] for i in self._keys]
    
    def __iter__(self):
        return self.keys()
        
    def __missing__(self, key):
        if key == 'LAT':
            j = arange(self._parent.dimensions['J'] + 1)
            data = j * self._parent.modelres[1] - 90
            if self._parent.halfpolar == 1:
                data[1:] -= 2
                data[-1] -= 2
            dims = ('J',)
            dtype = 'i'
            units = 'degrees north'
        elif key == 'LON':
            i = arange(self._parent.dimensions['I'] + 1)
            xres = self._parent.modelres[0]
            data = i * xres - (180 + xres / 2.) * self._parent.center180
            dims = ('I',)
            dtype = 'i'
            units = 'degrees east'
        elif key == 'tau0':
            tmp_key = self._keys[0]
            ord = self._tracer_ords[tmp_key]
            ids = self._ids(ord)
            data = self._headers[ids]['f10']
            dims = ('T',)
            dtype = 'i'
            units = 'hours since 0 GMT 1/1/1985'
        elif key == 'tau1':
            tmp_key = self._keys[0]
            ord = self._tracer_ords[tmp_key]
            ids = self._ids(ord)
            data = self._headers[ids]['f11']
            dims = ('T',)
            dtype = 'i'
            units = 'hours since 0 GMT 1/1/1985'
        else:
            ord = self._tracer_ords[key]
            ids = self._ids(ord)
            dims = ('T', 'K', 'J', 'I')
            dtype = 'f'
            base_units = self._headers[ids][0]['f9']
            scale = self._tracer_data[ord]['SCALE']
            carbon = self._tracer_data[ord]['C']
            units = self._tracer_data[ord]['UNIT']
            if carbon != 1 and 'C' not in units:
                warn("Scaling by carbon, but unit does not indicate carbon")
            data = self._data[ids] * scale * carbon
            
        return PseudoNetCDFVariable(self._parent, key, dtype, dims, units = units, long_name = key, var_desc = key, values = data)

    def _ids(self, ord):
        ids = []
        for hi, hd in enumerate(self._headers):
            if hd[8] == ord:
                ids.append(hi)
        
        if len(ids) == 0:
            raise KeyError
        else:
            return ids
            
    def _prune(self):
        ords = set([hd[8] for hd in self._headers])
        self._keys = [self._tracer_names[ord] for ord in ords]
        self._tracer_ords = dict([(self._tracer_names[ord], ord) for ord in ords])
                
    def __del__(self):
        del self._headers, self._data

    def __repr__(self):
        out = "{"
        for k in self._keys:
            out += "\n  '%s': PseudoNetCDFVariable(...)" % k
        out += "\n}"
        return out

    def __str__(self):
        return self.__repr__()


class bpch(PseudoNetCDFFile):
    """
    NetCDF-like class to interface with GEOS-Chem binary punch files
    
    f = bpch(path_to_binary_file, path_to_meta)
    dim = f.dimensions[dkey]
    var = f.variables[vkey]
    unit = var.unit
    dims = var.dimensions
    shape = var.shape
    """
    def __init__(self, bpch_path, tracerinfo = None, mode = 'r'):
        """
        bpch_path: path to binary punch file
        tracerinfo: path to ascii file with tracer definitions
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
        tracerinfo = tracerinfo or os.path.join(os.path.dirname(bpch_path), 'tracerinfo.dat')
        # Read binary data for general header and first datablock header
        header_block = fromfile(bpch_path, dtype = 'bool', count = _first_header_size)
        
        # combine data for convenience
        header = tuple(header_block[:_general_header_type.itemsize].view(_general_header_type)[0]) + \
                 tuple(header_block[_general_header_type.itemsize:].view(_datablock_header_type)[0])
        
        # Verify that all Fortran unformatted buffers match 
        try:
            assert(header[0] == header[2])
            assert(header[3] == header[5])
            assert(header[6] == header[11])
            assert(header[12] == header[-1])
        except AssertionError:
            raise ValueError, "BPCH Files fails header check"
        
        # Assign data from header to global attributes
        self.ftype = header[1]
        self.toptitle = header[4]
        self.modelname, self.modelres, self.halfpolar, self.center180 = header[7:11]
        self.category, self.tracer, self.unit, self.tau0, self.tau1, self.reserved, self.dim, self.start, self.skip = header[13:-1]
        
        # set datatype to appropriate shape
        data_type = dtype('>i4, %s>f4, >i4' % str(tuple(self.dim[::-1])))
        
        # set data block to block header and data
        block_type = dtype([('header', _datablock_header_type), ('data', data_type)])
        
        # load all data blocks
        self._memmap = memmap(bpch_path, dtype = block_type, offset = _general_header_type.itemsize, mode = mode)
        
        # Create variables and dimensions
        self.dimensions = dict(zip('I J K'.split(), self.dim))
        self.variables = _tracer_lookup(parent = self, headers = self._memmap['header'], data = self._memmap['data']['f1'], tracerinfo = tracerinfo)
        self.dimensions['T'] = self.variables['tau0'].shape[0]

    def __del__(self):
        del self.variables, self._memmap
      
    def __repr__(self):
        return PseudoNetCDFFile.__repr__(self) + str(self.variables)
if __name__ == '__main__':
    from numpy import median, indices, arange, meshgrid

    # Example: file open and variable aquisition
    f = bpch('restart.geos5.2005080100')
    var_key = 'Ox'
    var = f.variables[var_key]
    
    # Example: variable metadata print
    print var.long_name
    print var.dimensions
    print var.shape
    
    # Example: time/layer slice
    layer1 = var[0, 0, :, :]
    
    # Example: data statistics
    print layer1.min(), layer1.mean(), median(layer1), layer1.max(), var.units.strip()
    if False:
        # Example: spatial plotting
        from matplotlib import use
        use('Agg')
        from pylab import figure, xticks, yticks, title, colorbar
        from mpl_toolkits.basemap import Basemap
        
        # I'm actually not sure what the right projection for this 
        # data is.  So the next few lines might need to change.
        m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-180 * f.center180 - f.modelres[1] / 2.,\
                    urcrnrlon=360 - 180 * f.center180 - f.modelres[1] / 2.,\
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
        poly = m.pcolor(x, y, layer1)
        cb = colorbar(poly, ax = m.ax)
        cb.ax.set_xlabel(var.units.strip())
        xticks(meridians)
        yticks(parallels)
        title('Projection might be wrong')
        fig.savefig('layer1_%s.png' % var_key)
