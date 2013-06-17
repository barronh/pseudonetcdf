__all__ = ['bpch']
try:
    from bpch import bpch
except:
    import os

    # part of the default Python distribution
    from collections import defaultdict
    from warnings import warn

    # numpy is a very common installed library
    from numpy import fromfile, memmap, dtype, arange, zeros, ceil, diff, concatenate, append, pi, sin

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
            return defaultdict.__delitem__(self, key)
    
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
                raise KeyError("%s not found" % (key, ))

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
            self.variables = defaultpseudonetcdfvariable(list(groupvariables) + ['time', 'lev', 'tau0', 'tau1', 'crs', 'AREA', 'lat', 'lon', 'lat_bnds', 'lon_bnds'], getvar)
    
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
            self._special_keys = ['time', 'lev', 'tau0', 'tau1', 'crs', 'lat', 'lon', 'lat_bnds', 'lon_bnds']
            self._keys = keys + self._special_keys
            self._example_key = keys[0]
        
        def __missing__(self, key):
            if key in ('lat', 'lat_bnds'):
                yres = self._parent.modelres[1]
                if self._parent.halfpolar == 1:
                    data = concatenate([[-90.], arange(-90. + yres / 2., 90., yres), [90.]])
                else:
                    data = arange(-90, 90 + yres, yres)
            
                dims = ('lat',)
                dtype = 'i'
                kwds = dict(units = 'degrees north', long_name = key, var_desc = key)
                if key == 'lat':
                    data = data[:-1] + diff(data) / 2.
                    kwds['bounds'] = 'lat_bnds'
                else:
                    dims += ('nv',)
                    data = data.repeat(2,0)[1:-1].reshape(-1, 2)
                example = self[self._example_key]
                sj = getattr(example, 'STARTJ', 0)
                data = data[sj:sj + example.shape[2]]
                kwds = dict(standard_name = "latitude", long_name = "latitude", units = "degrees_north", axis = "Y", bounds = "lat_bnds")
            elif key in ('lon', 'lon_bnds'):
                xres = self._parent.modelres[0]
                i = arange(0, 360 + xres, xres)
                data = i - (180 + xres / 2. * self._parent.center180)
                dims = ('lon',)
                dtype = 'i'
                kwds = dict(units = 'degrees east', long_name = key, var_desc = key)
                if key == 'lon':
                    data = data[:-1] + diff(data) / 2.
                    kwds['bounds'] = 'lon_bnds'
                else:
                    dims += ('nv',)
                    data = data.repeat(2,0)[1:-1].reshape(-1, 2)
                example = self[self._example_key]
                si = getattr(example, 'STARTI', 0)
                data = data[si:si + example.shape[3]]
                kwds = dict(standard_name = "longitude", long_name = "longitude", units = "degrees_east", axis = "X", bounds = "lon_bnds")
            elif key == 'AREA':
               lon = self['lon']
               xres = self._parent.modelres[0]
               nlon = 360. / xres
               latb = self['lat_bnds']
               Re = self['crs'].semi_major_axis
               latb = append(latb[:, 0], latb[-1, 1])
               latr = pi / 180. * latb
               data = 2. * pi * Re * Re / (nlon) * ( sin( latr[1:] ) - sin( latr[:-1] ) )
               data = data[:, None].repeat(lon.size, 1)
               kwds = dict(units = 'm**2')
               dtype = 'i'
               dims = ('J', 'I')
            elif key == 'crs':
              dims = ()
              kwds = dict(grid_mapping_name = "latitude_longitude",
                          semi_major_axis = 6375000.0,
                          inverse_flattening = 0)
              dtype = 'i'
              data = zeros(1, dtype = dtype)
            elif key == 'time':
                tmp_key = self._keys[0]
                data = self._memmap[tmp_key]['header']['f10'] + .5
                dims = ('time',)
                dtype = 'i'
                kwds = dict(units = 'hours since 0 GMT 1/1/1985', standard_name = key, long_name = key, var_desc = key)
            elif key == 'lev':
                tmp_key = self._keys[0]
                data = arange(self._parent.dimensions['lev'], dtype = 'i')
                dims = ('lev',)
                dtype = 'i'
                kwds = dict(units = 'model layer', standard_name = key, long_name = key, var_desc = key)
            elif key == 'tau0':
                tmp_key = self._keys[0]
                data = self._memmap[tmp_key]['header']['f10']
                dims = ('time',)
                dtype = 'i'
                kwds = dict(units = 'hours since 0 GMT 1/1/1985', standard_name = key, long_name = key, var_desc = key)
            elif key == 'tau1':
                tmp_key = self._keys[0]
                data = self._memmap[tmp_key]['header']['f11']
                dims = ('time',)
                dtype = 'i'
                kwds = dict(units = 'hours since 0 GMT 1/1/1985', standard_name = key, long_name = key, var_desc = key)
            else:
                dtype = 'f'
                header = self._memmap[key]['header'][0]
                sl, sj, si = header['f14'][::-1] - 1
                group = header['f7'].strip()
                offset = self._diag_data[group]['offset']
                ord = header['f8'] + offset
                base_units = header['f9']
                scale = self._tracer_data[ord]['SCALE']
                carbon = self._tracer_data[ord]['C']
                units = self._tracer_data[ord]['UNIT']
                kwds = dict(scale = scale, carbon = carbon, units = units, base_units = base_units, standard_name = key, long_name = key, var_desc = key, coordinates = "time lev lat lon")
                tmp_data = self._memmap[key]['data']
                dims = ('time', 'lev', 'lat', 'lon')
                if 'srf_lev' in self._parent.dimensions:
                    if tmp_data.dtype['f1'].shape[0] == self._parent.dimensions['srf_lev']:
                        dims = ('time', 'srf_lev', 'lat', 'lon')
                
                assert((tmp_data['f0'] == tmp_data['f2']).all())
                data = tmp_data['f1'] * scale
                if any([sl != 0, sj != 0, si != 0]):
                    warn("%s is a subset variable" % key)
                    nl, nj, ni = header['f13'][::-1]
                    #import pdb; pdb.set_trace()
                    #tmp_data = zeros((data.shape[0], self._parent.dimensions['lev'], self._parent.dimensions['lat'], self._parent.dimensions['lon']), dtype = data.dtype)
                    #el, ej, ei = data.shape[1:]
                    #el += sl
                    #ej += sj
                    #ei += si
                    #tmp_data[:, sl:el, sj:ej, si:ei] = data[:]
                    #data = tmp_data
                    kwds['STARTI'] = si 
                    kwds['STARTJ'] = sj
                    kwds['STARTK'] = sl
            return PseudoNetCDFVariable(self._parent, key, dtype, dims, values = data, **kwds)
            
    class bpch(PseudoNetCDFFile):
        """
        NetCDF-like class to interface with GEOS-Chem binary punch files
    
        f = bpch(path_to_binary_file)
        dim = f.dimensions[dkey] # e.g., dkey = 'lon'
    
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

        def __init__(self, bpch_path, tracerinfo = None, diaginfo = None, mode = 'r', timeslice = slice(None)):
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
             timeslice: If the file is larger than 2GB, timeslice provides a way to subset results.
                        The subset requested depends on the data type of timeslice:
                            - int: return the a part of the file if it was broken into 2GB chunks (0..N-1)
                            - slice: return the times that correspond to that slice (i.e., range(ntimes)[timeslice])
                            - list/tuple/set: return specified times where each time is in the set (0..N-1)
            """
            self._ncattrs = () 
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
            self.dimensions = dict(zip('lon lat lev'.split(), dim))
            self.createDimension('nv', 2)
            tracerinfo = tracerinfo or os.path.join(os.path.dirname(bpch_path), 'tracerinfo.dat')
            if os.path.exists(tracerinfo):
                if os.path.isdir(tracerinfo): tracerinfo = os.path.join(tracerinfo, 'tracerinfo.dat')
                tracer_data = dict([(int(l[52:61].strip()), dict(NAME = l[:8].strip(), FULLNAME = l[9:39].strip(), MOLWT = float(l[39:49]), C = int(l[49:52]), TRACER = int(l[52:61]), SCALE = float(l[61:71]), UNIT = l[72:].strip())) for l in file(tracerinfo).readlines() if l[0] not in ('#', ' ')])
                tracer_names = dict([(k, v['NAME']) for k, v in tracer_data.iteritems()])
            else:
                warn('Reading file without tracerinfo.dat means that names and scaling are unknown')
                tracer_data = defaultdict(lambda: dict(SCALE = 1., C = 1.))
                tracer_names = defaultdictfromkey(lambda key: key)
        
            diaginfo = diaginfo or os.path.join(os.path.dirname(bpch_path), 'diaginfo.dat')
            if os.path.exists(diaginfo):
                if os.path.isdir(diaginfo): diaginfo = os.path.join(diaginfo, 'diaginfo.dat')
                diag_data = dict([(l[9:49].strip(), dict(offset = int(l[:8]), desc = l[50:].strip())) for l in file(diaginfo).read().strip().split('\n') if l[0] != '#'])
            else:
                warn('Reading file without diaginfo.dat loses descriptive information')
                diag_data = defaultdictfromkey(lambda key: dict(offset = 0, desc = key))
            
            if len(tracer_names) == 0 and not isinstance(tracer_names, defaultdictfromkey):
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
                tracer_number = header[8]
                unit = header[9].strip()
                if not isinstance(diag_data, defaultdictfromkey):
                    goffset = diag_data.get(group, {}).get('offset', 0)
                    try:
                        tracername = tracer_names[tracer_number + goffset]
                    except:
                        # There are some cases like with the adjoint where the tracer does not have
                        # a tracerinfo.dat line.  In this case, the name matches the tracer with that
                        # number (no offset).  The scaling, however, is not intended to be used.
                        # The unit for adjoint, for instance, is unitless.
                        if tracer_number not in tracer_names:
                            tracername = str(tracer_number)
                            tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 0, UNIT = unit)
                        else:
                            tracername = tracer_names[tracer_number]
                            tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = tracer_data[tracer_number]['C'], UNIT = unit)

                else:
                    warn('%s is not in diaginfo.dat; names and scaling cannot be resolved' % group)
                    goffset = 0
                    tracername = tracer_number
                    diag_data[group] = dict(offset = 0, desc = group)
                    tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 1., UNIT = unit)

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
            field_shapes = set([v[0].fields['data'][0].fields['f1'][0].shape for k, v in time_type.fields.iteritems()])
            field_levs = set([s_[0] for s_ in field_shapes])
            field_rows = set([s_[1] for s_ in field_shapes])
            field_cols = set([s_[2] for s_ in field_shapes])
            if len(field_levs) == 1:
                pass
            elif len(field_levs) == 2:
                self.dimensions['lev'] = max(field_levs)
                self.dimensions['srf_lev'] = min(field_levs)
            else:
                raise ValueError('Unclear how to handle more than 2 distinct layer sets')

            assert((float(os.path.getsize(bpch_path)) - _general_header_type.itemsize) % time_type.itemsize == 0.)
            # load all data blocks  
            try:
                datamap = memmap(bpch_path, dtype = time_type, offset = _general_header_type.itemsize, mode = mode)
            except OverflowError:
                hdrsize = _general_header_type.itemsize
                items = (2*1024**3-hdrsize) // time_type.itemsize
                if timeslice != slice(None):
                    filesize = os.stat(bpch_path).st_size
                    datasize = (filesize - hdrsize)
                    all_times = range(datasize / time_type.itemsize)
                    if isinstance(timeslice, int):
                        timeslice = slice(items * timeslice, items * (timeslice + 1))
                    if isinstance(timeslice, (list, tuple, set)):
                        times = timeslice
                    else:
                        times = all_times[timeslice]

                    outpath = bpch_path + '.tmp.part'
                    mint = times[0]
                    maxt = times[-1]
                    nt = maxt - mint + 1
                
                    if nt > items:
                        warn('Requested %d items; only returning %d items due to 2GB limitation' % (nt, items))
                        times = times[:items]

                    outfile = file(outpath, 'w')
                    infile = file(bpch_path, 'r')
                    hdr = infile.read(hdrsize)
                    outfile.write(hdr)
                    for t in all_times:
                        if t in times:
                            outfile.write(infile.read(time_type.itemsize))
                        else:
                            infile.seek(time_type.itemsize, 1)
                    outfile.close()
                    #print mint, maxt, nt, nt * time_type.itemsize
                    #cmd = 'dd if=%s ibs=%d skip=1 obs=%d | dd of=%s bs=%d skip=%d count=%d' % (bpch_path, hdrsize, time_type.itemsize, outpath, time_type.itemsize, mint, nt)
                    #print cmd
                    #os.system(cmd)
                    datamap = memmap(outpath, dtype = time_type, mode = mode, offset = hdrsize)
                else:
                    datamap = memmap(bpch_path, dtype = time_type, shape = (items,), offset = _general_header_type.itemsize, mode = mode)
                    warn('Returning only the first 2GB of data')

        
            # Create variables and dimensions
            self.variables = _tracer_lookup(parent = self, datamap = datamap, tracerinfo = tracer_data, diaginfo = diag_data, keys = keys)
            del datamap
            self.createDimension('time', self.variables['tau0'].shape[0])
            self.groups = dict([(k, _diag_group(self, k, v)) for k, v in self._groups.iteritems()])
            self.Conventions = "CF-1.6"


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
            path_to_test_file = raw_input('Enter path to a valid GEOS-Chem file:\n')

    group_key = f.groups.keys()[0]
    g = f.groups[group_key]
    i = 0; var_key = 'lat'
    while var_key in 'lat lon tau0 tau1'.split():
        var_key = g.variables.keys()[i]
        i += 1
    if 'Ox' in g.variables.keys():
        var_key = 'Ox'
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

        lat = f.variables['lat']
        lon = f.variables['lon']
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

