readers = []
def register(name, reader):
    global readers
    readers.append((name, reader))

def getreader(path):
    """
    Tests readers on path
    """
    for name, reader in readers:
        print name
        if hasattr(reader, 'isMine'):
            if reader.isMine(path):
                return reader
        else:
            try:
                reader(path)
                return reader
            except:
                pass
    else:
        raise IOError("Type was not identified")

def getfile(path):
    """
    Uses reader from getreader on path to return a file
    """
    reader = getreader(path)
    return reader(path)

from camxfiles.Memmaps import __all__ as camxmemmaps
from camxfiles.Memmaps import *
for reader in camxmemmaps:
    register(reader, eval(reader))

from geoschemfiles import __all__ as geoschem
from geoschemfiles import *
for reader in geoschem:
    register(reader, eval(reader))

from icarttfiles.ffi1001 import ffi1001
register('ffi1001', ffi1001)