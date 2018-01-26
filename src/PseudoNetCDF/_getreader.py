__all__ = ['getreader']

import unittest
import os
from warnings import warn
from PseudoNetCDF.netcdf import NetCDFFile

_readers = []
def testreader(reader, *args, **kwds):
    try:
        reader(*args, **kwds)
        return True
    except:
        return False

def getreader(*args, **kwds):
    """
    args - arguments for opening file
    kwds - keywords for file opener and optional format
    format - name of reader (optional)
    """
    global _readers
    format = kwds.pop('format', None)
    if not os.path.isfile(args[0]):
        warn('The first argument (%s) does not exist as a file.  First arguments are usually paths' % (args[0],))
    
    if format is None:
        _myreaders = _readers
        # Try to give preferential treatment based on suffix
        try:
            if len(args) > 0:
                ext = os.path.splitext(args[0])[1][1:]
            else:
                ext = os.path.splitext(kwds['path'])[1][1:]
            rdict = getreaderdict()
            if ext in rdict:
                _myreaders.insert(0, (ext, rdict[ext]))
        except:
            pass
    else:
        _myreaders = [(k, v) for k, v in _readers if format == k]
    
    for rn, reader in _myreaders:
        if getattr(reader, 'isMine', lambda *args, **kwds: testreader(reader, *args, **kwds))(*args, **kwds):
            return reader
    else:
        raise TypeError('No reader could open a file with these arguments %s %s' % (args, kwds))

def registerreader(name, reader):
    global _readers
    if not name in [k for k, v in _readers]:
        _readers.insert(0, (name, reader))
        return True
    else:
        return False

def pncopen(*args, **kwds):
    """
    Open any PNC supported format using args and kwds, which
    are format specific. format is not passed to the reader
    
    format = None, addcf = False, 
    args - arguments for opening file
    kwds - keywords for file opener and optional format and addcf
    format - name of reader (not passed to reader)
    addcf - boolean to add CF conventions (not passed to reader; default: True)
    diskless - boolean to add CF conventions (not passed to reader; default: False)
    """
    format = kwds.pop('format', None)
    addcf = kwds.pop('addcf', False)
    diskless = kwds.pop('diskless', False)
    if format is None:
        reader =  getreader(*args, format = format, **kwds)
    else:
        reader = getreaderdict()[format]
    
    outfile = reader(*args, **kwds)
    if addcf:
        try:
            from .conventions.ioapi import add_cf_from_ioapi
            from .core._wrapnc import WrapPNC
            from .core._functions import getvarpnc
            if not diskless:
                outfile = WrapPNC(outfile)
            else:
                outfile = getvarpnc(outfile, None)
            add_cf_from_ioapi(outfile)
        except Exception as e:
            warn(str(e))
            pass
    return outfile

anyfile = pncopen

def getreaderdict():
    return dict(_readers)
    
class TestPNCOPEN(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
    def testPNCOPENFMT(self):
        import PseudoNetCDF.testcase
        for format, path in PseudoNetCDF.testcase.self_described_paths.items():
            print('Test open with ', format, path)
            f = pncopen(path, format = format)

    def testPNCOPENNOFMT(self):
        import PseudoNetCDF.testcase
        for format, path in PseudoNetCDF.testcase.self_described_paths.items():
            print('Test open unspecified ', path)
            f = pncopen(path)

