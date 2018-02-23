from __future__ import print_function, unicode_literals
import sys
if sys.version_info.major == 2:
    bytes = lambda x, enc: str(x)
else:
    unicode = str

from warnings import warn
import re
import numpy as np
from collections import defaultdict, OrderedDict


from ._files import PseudoNetCDFFile
from ._variables import PseudoNetCDFMaskedVariable, PseudoNetCDFVariable
from ..userfuncs import *

import datetime

def pncrename(ifile, type_old_new):
    outf = getvarpnc(ifile, None)
    t,o,n = type_old_new.split(',')
    if t in ('d', 'dimension'):
        outf.dimensions[n] = outf.dimensions[o]
        del outf.dimensions[o]
        for k, v in outf.variables.items():
            if o in v.dimensions:
                v.dimensions = tuple([n if d == o else d for d in v.dimensions])
    if t in ('v', 'variable'):
        outf.variables[n] = outf.variables[o]
        del outf.variables[o]
    
    return outf
    
_translator = {'-': '_', '$': 'S', ' ': '_', '+': '_add_', '(': '', ')': ''}

def manglenames(f, translator = _translator):
    outf = getvarpnc(f, None)
    varkeys = list(outf.variables.items())
    for k, var in varkeys:
        nk = k
        for olds, news in _translator.items():
            nk = nk.replace(olds, news)
        if nk != k:
            outf.variables[nk] = outf.variables[k]
            del outf.variables[k]
    return outf

def removesingleton(f, rd, coordkeys = []):
    outf = PseudoNetCDFFile()
    for propkey in f.ncattrs():
        setattr(outf, propkey, getattr(f, propkey))
    for dk, d in f.dimensions.items():
        unlim = d.isunlimited()
        ni = len(d)
        if ni != 1 or dk != rd:
            tempd = outf.createDimension(dk, ni)
            if unlim:
                tempd.setunlimited(True)
    
    for vk, v in f.variables.items():
        dims = tuple([dk for dk in v.dimensions if dk in outf.dimensions])
        sdims = tuple([dk for dk in enumerate(v.dimensions) if dk[1] not in outf.dimensions])[::-1]
        propd = dict([(pk, getattr(v, pk)) for pk in v.ncattrs()])
        ov = outf.createVariable(vk, v.dtype.char, dims, **propd)
        outvals = v[...]
        for di, dk in sdims:
            outvals = outvals.take(0, axis = di)
        
        ov[...] = outvals[...]
    return outf
    
def getvarpnc(f, varkeys, coordkeys = [], copy = True):
    coordkeys = set(coordkeys)
    if varkeys is None:
        varkeys = list(f.variables.keys())
        varkeys = [k for k in varkeys if k not in coordkeys]
    else:
        newvarkeys = list(set(varkeys).intersection(f.variables.keys()))
        newvarkeys = [varkey for varkey in varkeys if varkey in newvarkeys]
        oldvarkeys = list(set(varkeys))
        oldvarkeys.sort()
        skipping = set(oldvarkeys).difference(newvarkeys)
        if len(skipping):
            warn('Skipping %s' % ', '.join(skipping))
        varkeys = newvarkeys

    if hasattr(f, 'copy'):
        outf = f.copy(props = True, dimensions = False, variables = False, data = False)
    else:
        from PseudoNetCDF.sci_var import WrapPNC
        outf = WrapPNC(f).copy(props = True, dimensions = False, variables = False, data = False)
    
    for varkey in varkeys:
        try:
            var = eval(varkey, None, f.variables)
        except:
            var = f.variables[varkey]
        for dimk, dimv in zip(var.dimensions, var.shape):
            if dimk not in outf.dimensions:
                newdimv = outf.createDimension(dimk, dimv)
                if f.dimensions[dimk].isunlimited():
                    newdimv.setunlimited(True)
                coordkeys.add(dimk)
                try:
                    tempv = f.variables[dimk]
                    if hasattr(tempv, 'bounds'):
                        coordkeys.add(tempv.bounds.strip())
                except (ValueError, KeyError, AttributeError) as e:
                    pass
        for coordk in coordkeys:
            if coordk in f.dimensions and coordk not in outf.dimensions:
                newdimv = outf.createDimension(coordk, len(f.dimensions[coordk]))
                if f.dimensions[coordk].isunlimited():
                    newdimv.setunlimited(True)
    
        propd = dict([(k, getattr(var, k)) for k in var.ncattrs()])
        if hasattr(var[...], 'fill_value') and 'fill_value' not in propd:
            propd['fill_value'] = var[...].fill_value
        
        if 'values' in propd:
            propd['pvalues'] = propd['values']
            del propd['values']
        if 'name' in propd:
            if not 'standard_name' in propd:
                propd['standard_name'] = propd['name']
            del propd['name']
        vals = var[...]
        if copy:
            vals = vals.copy()
        outf.createVariable(varkey, var.dtype.char, var.dimensions, values = vals, **propd)
    for coordkey in coordkeys:
        if coordkey in f.variables.keys():
            coordvar = f.variables[coordkey]
            propd = dict([(k, getattr(coordvar, k)) for k in coordvar.ncattrs()])
            outf.createVariable(coordkey, coordvar.dtype.char, coordvar.dimensions, values = coordvar[...], **propd)
            for dk in coordvar.dimensions:
                if dk not in outf.dimensions:
                    dv = outf.createDimension(dk, len(f.dimensions[dk]))
                    dv.setunlimited(f.dimensions[dk].isunlimited())
    
    return outf


def interpvars(f, weights, dimension, loginterp = []):
    """
    f - PseudoNetCDFFile
    weights - weights for new dimensions from old dimension dim(new, old)
    dimension - which dimensions will be reduced
    loginterp - iterable of keys to interp on log scale
    """
    outf = PseudoNetCDFFile()
    outf.dimensions = f.dimensions.copy()
    if hasattr(f, 'groups'):
        outf.groups = OrderedDict()
        for grpk, grpv in f.groups.items():
            outf.groups[grpk] = interpvars(grpv, weights, dimension)
    
    oldd = f.dimensions[dimension]
    didx, = [i for i, l in enumerate(weights.shape) if len(oldd) == l]
    
    newd = outf.createDimension(dimension, weights.shape[didx - 1])
    newd.setunlimited(oldd.isunlimited())
    for vark, oldvar in f.variables.items():
        if dimension in oldvar.dimensions:
            dimidx = list(oldvar.dimensions).index(dimension)
            if hasattr(oldvar, '_FillValue'):
                kwds = dict(fill_value = oldvar._FillValue, missing_value = oldvar._FillValue)
            else:
                kwds = dict()
            newvar = outf.createVariable(vark, oldvar.dtype.char, oldvar.dimensions, **kwds)
            for ak in oldvar.ncattrs():
                setattr(newvar, ak, getattr(oldvar, ak))
            if len(weights.shape) <= len(oldvar.dimensions):
                weightslice = (None,) * (dimidx) + (Ellipsis,) + (None,) * len(oldvar.dimensions[dimidx + 1:])
            else:
                weightslice = slice(None)        
            varslice = (slice(None,),) * dimidx + (None,)
            weightsv = weights[weightslice]
            oldvarv = oldvar[varslice]
            if not (weightsv.ndim - 1) == oldvar.ndim:
                warn('Wrong number of dimensions for %s' % (vark,))
            elif vark in loginterp:
                logv = np.ma.exp((weightsv * np.ma.log(oldvarv)).sum(dimidx + 1))
                newvar[:] = logv
            else:
                linv = (weightsv * oldvarv).sum(dimidx + 1)
                newvar[:] = linv
        else:
            outf.variables[vark] = oldvar
    return outf

def extract_from_file(f, lonlatfs, unique = False, gridded = None, method = 'nn', passthrough = True):
    from ..coordutil import getlonlatcoordstr
    lonlatcoordstr = ""
    for lonlatf in lonlatfs:
        lonlatcoordstr += getlonlatcoordstr(lonlatf)
    return extract_lonlat(f, lonlatcoordstr, unique = unique, gridded = gridded, method = method, passthrough = passthrough)
    
def extract_lonlat(f, lonlat, unique = False, gridded = None, method = 'nn', passthrough = True):
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    try:
        from StringIO import StringIO as BytesIO
    except ImportError:
        from io import BytesIO
    import os
    outf = PseudoNetCDFFile()
    outf.dimensions = f.dimensions.copy()
    if hasattr(f, 'groups'):
        outf.groups = {}
        for grpk, grpv in f.groups.items():
            outf.groups[grpk] = extract(grpv, lonlat)
    
    p2p = Pseudo2NetCDF()
    p2p.verbose = 0
    p2p.addGlobalProperties(f, outf)

    if gridded is None:
        gridded = ('longitude' in f.dimensions and 'latitude' in f.dimensions) or \
                  ('COL' in f.dimensions and 'ROW' in f.dimensions) or \
                  ('x' in f.dimensions and 'y' in f.dimensions)
    if isinstance(lonlat, (str, )):
        lonlat = [lonlat]
    lonlatin = lonlat
    lonlatout = []
    for ll in lonlat:
        if isinstance(ll, (str, )):
            try:
                if os.path.exists(ll):
                    ll = open(ll, 'r').read().strip()
            except Exception as e:
                warn('Windows machines may have uncessary warnings; ' + str(e))
                
            lonlatout.append(ll)
    lonlat = ('/'.join(lonlatout))
    try:
        lons, lats = np.genfromtxt(BytesIO(bytes(lonlat.replace('/', '\n'), 'ASCII')), delimiter = ',').T
    except Exception as e:
        print(str(e))
        raise e
    outf.lonlatcoords = lonlat
    if not method == 'll2ij':
        longitude = f.variables['longitude'][:]
        latitude = f.variables['latitude'][:]
        latlon1d = longitude.ndim == 1 and latitude.ndim == 1
    if method == 'll2ij':
        lonidxs, latidxs = f.ll2ij(lons, lats)
        def extractfunc(v, thiscoords):
            newslice = tuple([{'ROW': latidxs, 'COL': lonidxs, 'latitude': latidxs, 'longitude': lonidxs, 'points': latidxs, 'PERIM': latidxs}.get(d, slice(None)) for d in thiscoords])
            if newslice == ():
                return v
            else:
                return v[:][newslice]
    elif method == 'nn':
        if latlon1d and gridded:
            latitude = latitude[(slice(None), None, None)]
            longitude = longitude[(None, slice(None), None)]
        else:
            latitude = latitude[Ellipsis, None]
            longitude = longitude[Ellipsis, None]
    
        lonlatdims = latitude.ndim - 1
        londists = longitude - lons[(None,) * lonlatdims]
        latdists = latitude - lats[(None,) * lonlatdims]
        totaldists = ((latdists**2 + londists**2)**.5)
        if latlon1d and not gridded:
            latidxs, = lonidxs, = np.unravel_index(totaldists.reshape(-1, latdists.shape[-1]).argmin(0), totaldists.shape[:-1])
        else:
            latidxs, lonidxs = np.unravel_index(totaldists.reshape(-1, latdists.shape[-1]).argmin(0), totaldists.shape[:-1])
        def extractfunc(v, thiscoords):
            newslice = tuple([{'latitude': latidxs, 'longitude': lonidxs, 'points': latidxs, 'PERIM': latidxs}.get(d, slice(None)) for d in thiscoords])
            if newslice == ():
                return v
            else:
                return v[:][newslice]
    elif method == 'KDTree':
        if latlon1d and gridded:
            longitude, latitude = np.meshgrid(longitude, latitude)
        from scipy.spatial import KDTree
        tree = KDTree(np.ma.array([latitude.ravel(), longitude.ravel()]).T)
        dists, idxs = tree.query(np.ma.array([lats, lons]).T)
        if latlon1d and not gridded:
            latidxs, = lonidxs, = np.unravel_index(idxs, latitude.shape)
        else:
            latidxs, lonidxs = np.unravel_index(idxs, latitude.shape)
        def extractfunc(v, thiscoords):
            newslice = tuple([{'latitude': latidxs, 'longitude': lonidxs, 'points': latidxs, 'PERIM': latidxs}.get(d, slice(None)) for d in thiscoords])
            return v[newslice]
    elif method in ('linear', 'cubic'):
        from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator
        if method == 'cubic':
            interpclass = CloughTocher2DInterpolator
        else:
            interpclass = LinearNDInterpolator
        if latlon1d and gridded:
            longitude, latitude = np.meshgrid(longitude, latitude)
        points = np.array([longitude.ravel(), latitude.ravel()]).T
        def extractfunc(v, thiscoords):
            if not 'latitude' in thiscoords or not 'longitude' in thiscoords:
                return v
            newshape = [dl if d not in ('latitude', 'longitude') else -1 for di, (d, dl) in enumerate(zip(thiscoords, v.shape)) ]
            i1 = newshape.index(-1)
            if newshape.count(-1) > 1:
                i2 = newshape.index(-1, i1 + 1)
                assert(i1 == (i2 - 1))
                newshape.pop(i2)
            i2df = interpclass(points, np.rollaxis(v.reshape(*newshape), i1, 0))
            out = np.rollaxis(np.ma.array([i2df(lon, lat) for lat, lon in zip(lats, lons)]), 0, len(newshape))
            return out
        latidxs = extractfunc(latitude, ('latitude', 'longitude'))
    elif method in ('cubic', 'quintic'):
        from scipy.interpolate import interp2d
        if latlon1d and gridded:
            longitude, latitude = np.meshgrid(longitude, latitude)
        def extractfunc(v, thiscoords):
            i2df = interp2d(latitude, longitude, v, kind = method)
            return np.ma.array([i2df(lat, lon) for lat, lon in zip(lats, lons)])
        latidxs = extractfunc(latitude, '')
    else:
        raise ValueError('method must be: nn, KDTree')
    if unique:
        tmpx = OrderedDict()
        for lon, lat, lonlatstr in zip(lonidxs, latidxs, outf.lonlatcoords.split('/')):
            if (lon, lat) not in tmpx:
                tmpx[(lon, lat)] = lonlatstr
        
        lonidxs, latidxs = np.array(tmpx.keys()).T
        outf.lonlatcoords_orig = outf.lonlatcoords
        outf.lonlatcoords = '/'.join([tmpx[k] for k in zip(lonidxs, latidxs)])
    
    for k, v in f.variables.items():
        try:
            coords = v.coordinates.split()
        except:
            # special case for ioapi
            coords = [{'ROW': 'latitude', 'COL': 'longitude'}.get(tempdk, tempdk) for tempdk in v.dimensions]
        dims = v.dimensions
        outf.createDimension('points', len(latidxs))
        if passthrough or 'longitude' in coords or 'latitude' in coords:
            try:
                del outf.variables[k]
            except:
                pass
            newdims = []
            if len(dims) != len(coords):
                thiscoords = dims
            else:
                thiscoords = coords
            for d, c in zip(dims, thiscoords):
                if d not in ('longitude', 'latitude') and c not in ('longitude', 'latitude'):
                    newdims.append(d)
                else:
                    if 'points' not in newdims:
                        newdims.append('points')
                        
            
            newdims = tuple(newdims)
            newv = extractfunc(v, thiscoords)
            
            propd = dict([(ak, getattr(v, ak)) for ak in v.ncattrs()])
            nv = outf.createVariable(k, v.dtype.char, newdims, values = newv, **propd)
            setattr(nv, 'coordinates', getattr(v, 'coordinates', ' '.join(coords)))
            for di, dk in enumerate(newdims):
                if dk not in outf.dimensions:
                    outf.createDimension(dk, nv.shape[di])
    return outf

extract = extract_lonlat

def mask_vals(f, maskdef, metakeys = 'time layer level latitude longitude time_bounds latitude_bounds longitude_bounds ROW COL LAY TFLAG ETFLAG'.split()):
    mtype = maskdef.split(',')[0]
    mval = ','.join(maskdef.split(',')[1:])
    if mtype == 'where':
        maskexpr = 'np.ma.masked_where(mask, var[:].view(np.ndarray))'
        mask = eval(mval, None, f.variables)
    else:
        maskexpr = 'np.ma.masked_%s(var[:], %s)' % (mtype, mval)
    for varkey, var in f.variables.items():
        if varkey not in metakeys:
            try:
                vout = eval(maskexpr)
                f.variables[varkey] = PseudoNetCDFMaskedVariable(f, varkey, var.dtype.char, var.dimensions, values = vout, **dict([(pk, getattr(var, pk)) for pk in var.ncattrs()]))
            except Exception as e:
                warn('Cannot mask %s: %s' % (varkey, str(e)))
    return f
    
def slice_dim(f, slicedef, fuzzydim = True):
    """
    variables have dimensions (e.g., time, layer, lat, lon), which can be subset using 
        slice_dim(f, 'dim,start,stop,stride')
        
    e.g., slice_dim(f, 'layer,0,47,5') would sample every fifth layer starting at 0
    """
    inf = f

    historydef = "slice_dim(f, %s, fuzzydim = %s); " % (slicedef, fuzzydim)
    slicedef = slicedef.split(',')
    slicedef = [slicedef[0]] + list(map(eval, slicedef[1:]))
    if len(slicedef) == 2:
        slicedef.append(slicedef[-1] + 1)
    slicedef = (slicedef + [None,])[:4]
    dimkey, dmin, dmax, dstride = slicedef
    if dimkey not in inf.dimensions:
        warn('%s not in file' % dimkey)
        return inf

    unlimited = inf.dimensions[dimkey].isunlimited()
    if fuzzydim:
        partial_check = [key for key in inf.dimensions if dimkey == key[:len(dimkey)] and key[len(dimkey):].isdigit()]
        for dimk in partial_check:
            inf = slice_dim(inf, '%s,%s,%s,%s' % (dimk, dmin, dmax, dstride))
    
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    p2p = Pseudo2NetCDF(verbose = 0)
    outf = PseudoNetCDFFile()
    p2p.addDimensions(inf, outf)
    p2p.addGlobalProperties(inf, outf)
    
    for varkey in inf.variables.keys():
        var = inf.variables[varkey]
        if dimkey not in var.dimensions:
            p2p.addVariable(inf, outf, varkey)
        else:
            axis = list(var.dimensions).index(dimkey)
            vout = var[...].swapaxes(0, axis)[dmin:dmax:dstride].swapaxes(0, axis)
        
            newlen = vout.shape[axis]
            newdim = outf.createDimension(dimkey, newlen)
            newdim.setunlimited(unlimited)
            outf.variables[varkey] = vout
        
    history = getattr(outf, 'history', '')
    history += historydef
    setattr(outf, 'history', history)

    return outf
    
def _getfunc(a, func):
    """
    Get an approriate function that takes one optional keyword (axis)
    """
    if not hasattr(func, '__call__'):
        if hasattr(a, func):
            outfunc = getattr(a, func)
        elif isinstance(a, np.ma.MaskedArray):
            outfunc = lambda axis = None, keepdims = True: getattr(np.ma, func)(a, axis = axis, keepdims = keepdims)
        elif hasattr(np, func):
            outfunc = lambda axis = None, keepdims = True: getattr(np, func)(a, axis = axis, keepdims = keepdims)
        else:
            outfunc = lambda axis = None, keepdims = True: eval(func)(a, axis = axis, keepdims = keepdims)
    else:
        outfunc = lambda axis = None, keepdims = True: np.apply_along_axis(func1d = func, axis = axis, arr = a, keepdims = keepdims)
    return outfunc
    
def reduce_dim(f, reducedef, fuzzydim = True, metakeys = 'time layer level latitude longitude time_bounds latitude_bounds longitude_bounds ROW COL LAY TFLAG ETFLAG'.split()):
    """
    variable dimensions can be reduced using
    
    reduce_dim(file 'dim,function,weight')
    
    e.g., reduce_dim(layer,mean,weight).
    
    Weighting is not fully functional.
    """
    inf = f
    metakeys = [k for k in metakeys if k in inf.variables.keys()]
    historydef = "reduce_dim(f, %s, fuzzydim = %s, metakeys = %s); " % (reducedef, fuzzydim, metakeys)
    import numpy as np
    if hasattr(reducedef, 'split') and hasattr(reducedef, 'count'):
        commacount = reducedef.count(',')
        reducevals = reducedef.split(',')
    else:
        commacount = len(reducedef)
        reducevals = reducedef
    if commacount == 3:
        dimkey, func, numweightkey, denweightkey = reducevals
        numweight = inf.variables[numweightkey]
        denweight = inf.variables[denweightkey]
    elif commacount == 2:
        dimkey, func, numweightkey = reducevals
        numweight = inf.variables[numweightkey]
        denweightkey = None
    elif commacount == 1:
        dimkey, func = reducevals
        numweightkey = None
        denweightkey = None
    if fuzzydim:
        partial_check = [key for key in inf.dimensions if dimkey == key[:len(dimkey)] and key[len(dimkey):].isdigit()]
        for dimk in partial_check:
            if commacount == 1:
                inf = reduce_dim(inf, '%s,%s' % (dimk, func),)
            elif commacount == 2:
                inf = reduce_dim(inf, '%s,%s,%s' % (dimk, func, numweightkey),)
            elif commacount == 3:
                inf = reduce_dim(inf, '%s,%s,%s,%s' % (dimk, func, numweightkey, denweightkey),)
    if dimkey not in inf.dimensions:
        warn('%s not in file' % dimkey)
        return inf

    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    p2p = Pseudo2NetCDF(verbose = 0)
    outf = PseudoNetCDFFile()
    p2p.addDimensions(inf, outf)
    del outf.dimensions[dimkey]
    p2p.addGlobalProperties(inf, outf)
    
    #unlimited = inf.dimensions[dimkey].isunlimited()
    #outf.createDimension(dimkey, 1)
    #if unlimited:
    #    outf.dimensions[dimkey].setunlimited(True)

    for varkey in inf.variables.keys():
        var = inf.variables[varkey]
        if dimkey not in var.dimensions:
            p2p.addVariable(inf, outf, varkey)
            continue
        
        axis = list(var.dimensions).index(dimkey)
        #def addunitydim(var):
        #    return var[(slice(None),) * (axis + 1) + (None,)]
        vreshape = var[slice(None)]
        #vreshape = addunitydim(var)
        if not varkey in metakeys:
            if numweightkey is None:
                vout = _getfunc(vreshape, func)(axis = axis, keepdims = True)
            elif denweightkey is None:
                wvar = var * np.array(numweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)]
                vout = getattr(wvar[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis)
                vout.units = vout.units.strip() + ' * ' + numweight.units.strip()
                if hasattr(vout, 'base_units'):
                    vout.base_units = vout.base_units.strip() + ' * ' + numweight.base_units.strip()
            else:
                nwvar = var * np.array(numweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]),)]
                vout = getattr(nwvar[(slice(None),) * (axis + 1) + (None,)], func)(axis = axis) / getattr(np.array(denweight, ndmin = var.ndim)[(slice(None),)*axis + (slice(0,var.shape[axis]), None)], func)(axis = axis)
        else:
            if '_bounds' not in varkey and '_bnds' not in varkey:
                vout = _getfunc(vreshape, func)(axis = axis, keepdims = True)
            else:
                vout = _getfunc(vreshape, func)(axis = axis, keepdims = True)
                vmin = _getfunc(vreshape, 'min')(axis = axis, keepdims = True)
                vmax = _getfunc(vreshape, 'max')(axis = axis, keepdims = True)
                if 'lon' in varkey or 'time' in varkey:
                    try:
                        vout[..., [0, 3]] = vmin[..., [0, 3]]
                        vout[..., [1, 2]] = vmax[..., [1, 2]]
                    except:
                        vout[..., [0, 1]] = vmin[0, 0], vmax[0, 1]
                elif 'lat' in varkey:
                    nmin = vout.shape[-1] // 2
                    vout[..., :nmin] = vmin[..., :nmin]
                    vout[..., nmin:] = vmax[..., nmin:]
        if dimkey not in outf.dimensions:
             outdim = outf.createDimension(dimkey, vout.shape[axis])
             outdim.setunlimited(inf.dimensions[dimkey].isunlimited())
        nvar = outf.variables[varkey] = PseudoNetCDFMaskedVariable(outf, varkey, var.dtype.char, var.dimensions, values = vout)
        for k in var.ncattrs():
            setattr(nvar, k, getattr(var, k))

    history = getattr(outf, 'history', '')
    history += historydef
    setattr(outf, 'history', history)
    return outf

def pncfunc(func, ifile1, coordkeys = [], verbose = 0):
    """
    Perform function (func) on all variables in ifile1.  The returned file (rfile) contains the result
    
    rfile = ifile1 <op>
    
    func can be a function or string
    """
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    
    # Copy infile1 to a temporary PseudoNetCDFFile
    p2p = Pseudo2NetCDF()
    p2p.verbose = verbose
    tmpfile = PseudoNetCDFFile()
    p2p.convert(ifile1, tmpfile)
    
    # For each variable, assign the new value
    # to the tmpfile variables.
    for k in tmpfile.variables.keys():
        if k in coordkeys: continue
        outvar = tmpfile.variables[k]
        in1var = ifile1.variables[k]
        if not hasattr(func, '__call__'):
            if hasattr(in1var, func):
                outval = getattr(in1var, func)()
            elif '.' == func[:1]:
                outval = eval('in1var[:]' + func)
        else:
            outval = func(in1var[:])
        outval = np.ma.filled(np.ma.masked_invalid(outval), -999)
        if outvar.ndim > 0:
            outvar[:] = outval
        else:
            outvar.itemset(outval)
        outvar.fill_value = -999
    return tmpfile

def pncbo(op, ifile1, ifile2, coordkeys = [], verbose = 0):
    """
    Perform binary operation (op) on all variables in ifile1
    and ifile2.  The returned file (rfile) contains the result
    
    rfile = ifile1 <op> ifile2
    
    op can be any valid operator (e.g., +, -, /, *, **, &, ||)
    """
    tmpfile = ifile1.copy(props = True, dimensions = True, variables = False, data = False)
    
    # For each variable, assign the new value
    # to the tmpfile variables.
    for k in ifile1.variables.keys():
        in1var = ifile1.variables[k]
        if k not in ifile2.variables.keys() or k in coordkeys:
            warn('%s not found in ifile2' % k)
            tmpfile.copyVariable(in1var, newkey = k)
        else:
            in2var = ifile2.variables[k]
            propd = dict([(ak, getattr(in1var, ak)) for ak in in1var.ncattrs()])
            unit1 = getattr(in1var, 'units', 'unknown')
            unit2 = getattr(in2var, 'units', 'unknown')
            propd['units'] = '(%s) %s (%s)' % (unit1, op, unit2)
            outval = np.ma.masked_invalid(eval('in1var[...] %s in2var[...]' % op).view(np.ndarray))
            outvar = tmpfile.createVariable(k, in1var.dtype.char, in1var.dimensions, fill_value = -999, values = outval)
            outvar.setncatts(propd)
    return tmpfile

def pncbfunc(func, ifile1, ifile2, coordkeys = [], verbose = 0):
    """
    Perform binary function (func) on all variables in ifile1
    and ifile2.  The returned file (rfile) contains the result
    
    rfile = ifile1 <op> ifile2
    
    op can be any valid operator (e.g., +, -, /, *, **, &, ||)
    """
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    
    # Copy infile1 to a temporary PseudoNetCDFFile
    p2p = Pseudo2NetCDF()
    p2p.verbose = verbose
    tmpfile = PseudoNetCDFFile()
    p2p.convert(ifile1, tmpfile)
    
    # For each variable, assign the new value
    # to the tmpfile variables.
    for k in tmpfile.variables.keys():
        if k in coordkeys: continue
        outvar = tmpfile.variables[k]
        in1var = ifile1.variables[k]
        if k not in ifile2.variables.keys():
            warn('%s not found in ifile2' % k)
            continue
        in2var = ifile2.variables[k]
        outval = np.ma.filled(np.ma.masked_invalid(func(in1var[...], in2var[...])), -999)
        if outvar.ndim > 0:
            outvar[:] = outval
        else:
            outvar.itemset(outval)
        outvar.fill_value = -999
    return tmpfile

def _namemangler(k):
    k = k.replace('$', 'dollar')
    k = k.replace('-', 'hyphen')
    k = k.replace('(', 'lparen')
    k = k.replace(')', 'rparen')
    return k

def pncexpr(expr, ifile, verbose = 0):
    """
    Evaluate an arbitrary expression in the context of ifile.variables
    and add the result to the file with appropriate units.
    """
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    from PseudoNetCDF.sci_var import WrapPNC
    from symtable import symtable
    
    # Copy file to temporary PseudoNetCDF file
    comp = compile(expr, 'none', 'exec')

                    
    #varkeys = [key for key in comp.co_names if key in ifile.variables]
    #getvarpnc(ifile, varkeys)
    varpnc = tmpfile = WrapPNC(ifile, None)

    # Get NetCDF variables as a dictionary with 
    # names mangled to allow special characters
    # in the names
    vardict = dict([(_namemangler(k), varpnc.variables[k]) for k in varpnc.variables.keys()])
    # Compile the expr
    symtbl = symtable(expr, '<pncexpr>', 'exec')
    symbols = symtbl.get_symbols()
    for symbol in symbols:
        key = symbol.get_name()
        if key in vardict:
            tmpvar = vardict[key]
            break
    else:
        tmpvar = PseudoNetCDFVariable(None, 'temp', 'f', ())
    
    propd = dict([(k, getattr(tmpvar, k)) for k in tmpvar.ncattrs()])
    dimt = tmpvar.dimensions
    # Add all used constants as properties
    # of the output file
    vardict['ifile'] = ifile
    vardict['infile'] = ifile
    vardict['np'] = np
    exec('from scipy.constants import *', None, vardict)
    for k in ifile.ncattrs():
        if not k in vardict:
            vardict[k] = getattr(ifile, k)
    oldkeys = set(vardict.keys())
    
    # Assign expression to new variable.
    exec(comp, None, vardict)
    
    #newkeys = set(vardict.keys())
    #assignedkeys = [k for k in newkeys.difference(oldkeys) if k[:1] != '_']
    assignedkeys = [s.get_name() for s in symbols if s.is_assigned()]
    assignedkeys = [k for k in assignedkeys if k in vardict]
    for key in assignedkeys:
        val = vardict[key]
        # if the output variable has no dimensions, there is likely a problem
        # and the output should be defined.
        if isinstance(val, (PseudoNetCDFVariable,)) and val.dimensions != ():
            tmpfile.variables[key] = val
        else:
            tmpfile.createVariable(key, val.dtype.char, dimt, values = val, **propd)
    
    return tmpfile
    
def seqpncbo(ops, ifiles, coordkeys = []):
    for op in ops:
        ifile1, ifile2 = ifiles[:2]
        newfile = pncbo(op = op, ifile1 = ifile1, ifile2 = ifile2, coordkeys = coordkeys)
        del ifiles[:2]
        ifiles.insert(0, newfile)
    return ifiles

def mesh_dim(f, mesh_def):
    dimkey, meshfactor, aggfunc = mesh_def.split(',')
    meshfactor = float(meshfactor)
    spread=lambda a, n, axis: a.repeat(n, axis) * meshfactor
    try:
        aggfunc = eval(aggfunc)
    except:
        aggfunc = getattr(np, aggfunc)
    if meshfactor < 1.:
        oldres = int(1./meshfactor)
        assert(1./meshfactor == oldres)
        newres = 1
    elif meshfactor > 1.:
        newres = int(meshfactor)
        assert(meshfactor == newres)
        oldres = 1
    from PseudoNetCDF.MetaNetCDF import newresolution
    nrf = newresolution(f, dimension = dimkey, oldres = oldres, newres = newres, repeat_method = aggfunc, condense_method = aggfunc)
    f.dimensions[dimkey] = nrf.dimensions[dimkey]
    for k, v in f.variables.items():
        if dimkey in v.dimensions:
            f.variables[k] = nrf.variables[k]
    return f

def add_attr(f, attr_def):
    pieces = attr_def.split(',')
    att_nm, var_nm, mode, att_typ = pieces[:4]
    att_val = ','.join(pieces[4:])
    if var_nm == 'global':
        var = f
    else:
        var = f.variables[var_nm]

    if not (att_typ == 'c' and isinstance(att_val, (str, unicode))):
        att_val = np.array(eval(att_val), dtype = att_typ)

    if mode in ('a',):
        att_val = np.append(getattr(var, att_nm, []), att_val)

    if mode in ('a', 'c', 'm', 'o'):
        setattr(var, att_nm, att_val)
    elif mode in ('d',):
        delattr(var, att_nm)
    else:
        raise KeyError('mode must be either a c m o or d')

def convolve_dim(f, convolve_def):
    convolve_parts = convolve_def.split(',')
    dimkey = convolve_parts.pop(0)
    mode = convolve_parts.pop(0)
    weights = np.array(convolve_parts, dtype = 'f')
    outf = PseudoNetCDFFile()
    from PseudoNetCDF.pncgen import Pseudo2NetCDF
    p2p = Pseudo2NetCDF(verbose = 0)
    p2p.addGlobalProperties(f, outf)
    p2p.addDimensions(f, outf)
    dim = outf.dimensions[dimkey]
    dim = outf.createDimension(dimkey, len(np.convolve(weights, np.arange(len(dim)), mode = mode)))
    dim.setunlimited(f.dimensions[dimkey].isunlimited())
    for vark, var in f.variables.items():
        lconvolve = dimkey in var.dimensions
        p2p.addVariable(f, outf, vark, data = not lconvolve)
        if lconvolve:
            axisi = list(var.dimensions).index(dimkey)
            values = np.apply_along_axis(func1d = lambda x_: np.convolve(weights, x_, mode = mode), axis = axisi, arr = var[:])
            if isinstance(var[:], np.ma.MaskedArray):
                values = np.ma.masked_invalid(values)
            
            outf.variables[vark][:] = values
    return outf

def merge(fs):
    outf = getvarpnc(fs[0], None)
    for f in fs[1:]:
        for p in f.ncattrs():
            if not p in outf.ncattrs():
                setattr(outf, p, getattr(f, p))
        for d, v in f.dimensions.items():
            if not d in outf.dimensions:
                nv = outf.createDimension(d, len(v))
                nv.setunlimited(v.isunlimited())
        for k, v in f.variables.items():
            if k in outf.variables:
                if v.shape != outf.variables[k].shape or not (v[...] == outf.variables[k][...]).all():
                    warn('%s already in output' % k)
            else:
                propd = dict([(p, getattr(v, p)) for p in v.ncattrs()])
                var = outf.createVariable(k, v.dtype.char, v.dimensions, values = v, **propd)
    
    return outf

def stack_files(fs, stackdim, coordkeys = []):
    """
    Create files with dimensions extended by stacking.
    
    Currently, there is no sanity check...
    
    """
    f = PseudoNetCDFFile()
    tmpf = fs[0]
    dimensions = [f_.dimensions for f_ in fs]
    shareddims = {}
    for dimk, dim in tmpf.dimensions.items():
        if dimk == stackdim:
            continue
        dimlens = map(len, [dims[dimk] for dims in dimensions])
        if all([len(dim) == i for i in dimlens]):
            shareddims[dimk] = len(dim)
    differentdims = [set(dims.keys()).difference(shareddims.keys()) for dims in dimensions]
    assert(all([different == set([stackdim]) for different in differentdims]))
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    p2p = Pseudo2NetCDF(verbose = 0)
    p2p.addDimensions(tmpf, f)
    f.createDimension(stackdim, sum([len(dims[stackdim]) for dims in dimensions]))
    p2p.addGlobalProperties(tmpf, f)
    for tmpf in fs:
        for varkey, var in tmpf.variables.items():
            if not stackdim in var.dimensions:
                if varkey in f.variables:
                    if not varkey in coordkeys:
                        warn('Got duplicate variables for %s without stackable dimension; first value retained' % varkey)
                else:
                    p2p.addVariable(tmpf, f, varkey, data = True)
            else:
                if not varkey in f.variables.keys():
                    axisi = list(var.dimensions).index(stackdim)
                    values = np.ma.concatenate([f_.variables[varkey][:] for f_ in fs], axis = axisi)
                    p2p.addVariable(tmpf, f, varkey, data = False)
                    f.variables[varkey][:] = values
        
    return f

def splitdim(inf, olddim, newdims, newshape):
    oldsize = len(inf.dimensions[olddim])
    newsize = np.prod(newshape)
    if newsize != oldsize:
        raise ValueError('New shape, must match old dimension length: %d %d %s' % (oldsize, newsize, newshape))
    if len(newdims) != len(newshape):
        raise ValueError('Shape and dimensions must match in length')
    from PseudoNetCDF.sci_var import Pseudo2NetCDF
    p2n = Pseudo2NetCDF()
    outf = PseudoNetCDFFile()
    for dk, d in inf.dimensions.items():
        if dk == olddim:
            for dk, dl in zip(newdims, newshape):
                outf.createDimension(dk, dl)
        else:
            p2n.addDimension(inf, outf, dk)
    
    for vk, invar in inf.variables.items():
        if olddim in invar.dimensions:
            outdims = []
            outshape = []
            for dk in invar.dimensions:
                if dk == olddim:
                    outdims.extend(newdims)
                    outshape.extend(newshape)
                else:
                    outdims.append(dk)
                    outshape.append(len(inf.dimensions[dk]))
            
            outvar = outf.createVariable(vk, invar.dtype.char, tuple(outdims))
            p2n.addVariableProperties(invar,outvar)
            outvar[:] = invar[:].reshape(*outshape)
        else:
            p2n.addVariable(inf, outf, vk)
    
    return outf
