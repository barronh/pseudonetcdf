from tempfile import NamedTemporaryFile as tnf

from _files import PseudoNetCDFFile
from PseudoNetCDF.netcdf import NetCDFFile

def get_ncf_object(path_or_object, mode, format = 'NETCDF4_CLASSIC'):
    """
    Return an open NetCDF or PseudoNetCDF from an object or path
    """
    from os.path import exists, isfile, isdir
    from PseudoNetCDF.netcdf import NetCDFFile
    read_only = ('r', 'r+', 'rs', 'rs+', 'r+s')
    if isinstance(path_or_object, (str, unicode)):
        if exists(path_or_object):
            if isfile(path_or_object):
                ncf_object = NetCDFFile(path_or_object, mode, format = format)
            elif isdir(path_or_object):
                raise ValueError("Got directory at %s; not sure what to do" % path_or_object)
            else:
                raise ValueError("Expected file or directory at %s" % path_or_object)
        elif mode not in read_only:
            ncf_object = NetCDFFile(path_or_object, mode, format = format)
        else:
            raise IOError("Cannot open missing file for reading")
    elif isinstance(path_or_object, NetCDFFile) or isinstance(path_or_object, PseudoNetCDFFile):
        return path_or_object
    elif path_or_object is None and mode not in read_only:
        tfile=tnf(mode='w+b')
        npath=tfile.name
        ncf_object=NetCDFFile(npath,mode)
    else:
        raise ValueError("Not a path; not a netCDF file; not a PseudoNetCDF file... I don't know what to do")
    return ncf_object

def get_dimension_length(pfile, key):
    """
    Return the length of a dimension (key) from pfile
    """
    dim = pfile.dimsensions[key]
    if dim is None:
        for k in pfile.variables.keys():
            v = pfile.variables[k]
            if key in v.dimensions:
                return v[...].shape[list(v.dimensions).index(key)]
        return 0
    elif isinstance(dim, int):
        return dim
    else:
        return len(dim)
    
