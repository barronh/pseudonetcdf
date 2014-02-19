__all__ = ['getreader']

import os
from warnings import warn
from PseudoNetCDF.netcdf import NetCDFFile

_readers = [('netcdf', NetCDFFile)]
def testreader(reader, *args, **kwds):
    try:
        reader(*args, **kwds)
        return True
    except:
        return False

def getreader(*args, **kwds):
    global _readers
    if not os.path.isfile(args[0]):
        warn('The first argument (%s) does not exist as a file.  First arguments are usually paths' % (args[0],))
    
    for rn, reader in _readers:
        if getattr(reader, 'isMine', lambda *args, **kwds: testreader(reader, *args, **kwds))(*args, **kwds):
            return reader
    else:
        raise TypeError('No reader could open a file with these arguments %s %s' % (args, kwds))

def registerreader(name, reader):
    global _readers
    _readers.insert(0, (name, reader))

def anyfile(*args, **kwds):
    return getreader(*args, **kwds)(*args, **kwds)

def getreaderdict():
    return dict(_readers)