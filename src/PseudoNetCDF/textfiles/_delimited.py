from PseudoNetCDF.sci_var import PseudoNetCDFFile
import numpy as np

class csv(PseudoNetCDFFile):
    def __init__(self, path, coordkeys = "time time_bounds TFLAG ETFLAG latitude latitude_bounds longitude longitude_bounds lat lat_bnds lon lon_bnds etam_pressure etai_pressure layer_bounds layer47 layer".split(), delimiter = ',', names = True, **kwds):
        """
        path - place to find csv file
        coordkeys - use these keys as dimensions and coordinate variables
        delimiter - use this as delimiter (default = ',')
        names - see help in recfromtxt (Default = True)
        kwds - np.recfromtxt keywords
        
        * Note: currently only works when all coordinate variables are 1-d
        """
        kwds['names'] = names
        kwds['delimiter'] =  delimiter
        data = np.recfromtxt(path, **kwds)
        dimkeys = [dk for dk in data.dtype.names if dk in coordkeys]
        varkeys = [vk for vk in data.dtype.names if not vk in coordkeys]
        for dk in dimkeys:
            dv = np.unique(data[dk])
            dv.sort()
            self.createDimension(dk, len(dv))
            dvar = self.createVariable(dk, dv.dtype.char, (dk,))
            dvar[:] = dv
        
        for vk in varkeys:
            vv = data[vk]
            var = self.createVariable(vk, vv.dtype.char, tuple(dimkeys))
            for idx in np.ndindex(var.shape):
                thisidx = np.sum([data[dk] == self.variables[dk][di] for di, dk in zip(idx, dimkeys)], axis = 0) == len(dimkeys)
                if thisidx.any():
                    var[idx] = vv[thisidx]
        
    

def ncf2csv(ifile, outpath, delimiter = ',', coordkeys = "time time_bounds TFLAG ETFLAG latitude latitude_bounds longitude longitude_bounds lat lat_bnds lon lon_bnds etam_pressure etai_pressure layer_bounds layer47 layer".split()):
    header = [k for k, v in ifile.variables.iteritems() if k not in coordkeys and v.size > 1 and k not in ifile.dimensions]
    dims = set([ifile.variables[k].dimensions for k in header])
    if len(dims) > 1:
        warn('Making multiple csv outputs because not all output variables have the same dimensions')
    dimsets = {}
    for di, dim in enumerate(dims):
        if len(dims) > 1:
            outfile = open(outpath + str(di), 'w')
        else:
            outfile = open(outpath, 'w')
        
        dimdict = [(di, k) for di, k in enumerate(dim) if k in ifile.variables]
        dimheader = [k for k in dim if k in ifile.variables]
        header = dimsets[dim] = [k for k, v in ifile.variables.iteritems() if v.dimensions == dim]
        dimvars = [ifile.variables[k] for k in dimheader]
        vars = [ifile.variables[k] for k in header]
        print >> outfile, delimiter.join(dimheader + header)
        for idx in np.ndindex(ifile.variables[header[-1]].shape):
            outvals = []
            for dk, dv in zip(dimheader, dimvars):
                dv = ifile.variables[dk]
                didx = tuple([iidx for i, iidx in enumerate(idx) if dim[i] in dv.dimensions])
                outvals.append(repr(dv[didx]))
            for vk, vv in zip(header, vars):
                outvals.append(repr(vv[idx]))
            
            print >> outfile, delimiter.join(outvals)
    
