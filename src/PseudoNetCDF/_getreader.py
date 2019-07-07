__all__ = ['getreader', 'pncopen', 'pncmfopen']

import unittest
import os
from warnings import warn
from PseudoNetCDF.netcdf import NetCDFFile

_readers = [('Dataset', NetCDFFile)]


def testreader(reader, *args, **kwds):
    try:
        reader(*args, **kwds)
        return True
    except Exception:
        return False


def getreader(*args, **kwds):
    """
    Parameters
    ----------
    *args :
        arguments for opening file
    **kwds :
        keywords for file opener and optional format
    format : string
        name of reader (optional)

    Returns
    -------
    reader : class
    """
    global _readers
    format = kwds.pop('format', None)
    if not os.path.isfile(args[0]):
        warn(('The first argument (%s) does not exist as a file.  ' +
              'First arguments are usually paths') % (args[0],))

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
        except Exception:
            pass
    else:
        _myreaders = [(k, v) for k, v in _readers if format == k]

    for rn, reader in _myreaders:
        if hasattr(reader, 'isMine'):
            checker = reader.isMine
        else:
            def checker(*args, **kwds):
                try:
                    testreader(reader, *args, **kwds)
                    return True
                except Exception:
                    return False
        if checker(*args, **kwds):
            return reader
    else:
        raise TypeError(('No reader could open a file with these ' +
                         'arguments %s %s') % (args, kwds))


def registerreader(name, reader):
    global _readers, pncopen
    if name not in [k for k, v in _readers]:
        _readers.insert(0, (name, reader))
        if pncopen.__doc__ is not None:
            newdoc = '\n        - '.join([pncopen.__doc__, name])
            pncopen.__doc__ = newdoc
        return True
    else:
        return False


def pncmfopen(*args, stackdim=None, **kwds):
    """Open any PNC supported format using pncopen on all files
    passed as teh first argument of pncmfopen. See pncopen

    Parameters
    ----------
    args : arguments
        args[0] must be a list of paths, other arguments are format specific
    stackdim : str
        dimension upon which to stack files
    kwds : dict
        see pncopen for more details

    Returns
    -------
    pfile : PseudoNetCDF
    """
    paths = args[0]
    files = [pncopen(path, *args[1:], **kwds) for path in paths]
    file1 = files[0]
    return file1.stack(files[1:], stackdim=stackdim)


def pncopen(*args, **kwds):
    """Open any PNC supported format using args and kwds, which
    are format specific. format is not passed to the reader

    Parameters
    ----------
    *args :
        arguments for opening file
    **kwds :
        keywords for reader (see format)
    format : string
        name of reader (not passed to reader), default auto-detect
        see Format Options for formats
    addcf : boolean
        to add CF conventions (not passed to reader; default: False)
    diskless : boolean
        to add CF conventions (not passed to reader; default: False)
    help : boolean
        without format, returns help of pncopen and with format keyword,
        returns help of that class. See the __init__ interface for help
        with keywords that are not in the class __doc__.

    Returns
    -------
    pfile : PseudoNetCDF

    Notes
    -----
    Format Options:
        - add your own by subclassing PseudoNetCDFFile
        - for a full list, use pncopen(help=True)
    """
    formathelp = kwds.pop('help', False)
    format = kwds.pop('format', None)
    addcf = kwds.pop('addcf', False)
    diskless = kwds.pop('diskless', False)
    if format is None:
        if formathelp:
            return help(pncopen)
        reader = getreader(*args, format=format, **kwds)
    else:
        reader = getreaderdict()[format]
        if formathelp:
            return help(reader)

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
            pncopen(path, format=format)

    def testPNCOPENNOFMT(self):
        import PseudoNetCDF.testcase
        for format, path in PseudoNetCDF.testcase.self_described_paths.items():
            print('Test open unspecified ', path)
            pncopen(path)
