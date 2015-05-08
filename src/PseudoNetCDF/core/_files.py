import unittest
from PseudoNetCDF._getreader import registerreader
from PseudoNetCDF.netcdf import NetCDFFile
from collections import OrderedDict
from _dimensions import PseudoNetCDFDimension
from _variables import PseudoNetCDFVariable, PseudoNetCDFMaskedVariable

class OrderedDefaultDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

class classreg(type):
    def __init__(cls, name, bases, clsdict):
        pieces = str(cls).split('\'')[1].split('.')
        longname = '.'.join([p for p in pieces[1:-1]  if '_' != p[0] and p not in ('core',)] + [pieces[-1]])
        if len(cls.mro()) > 2:
            if name != 'PseudoNetCDFFileMemmap':
                registerreader(name, cls)
                registerreader(longname, cls)
        super(classreg, cls).__init__(name, bases, clsdict)

class PseudoNetCDFFile(object):
    """
    PseudoNetCDFFile provides an interface and standard set of
    methods that a file should present to act like a netCDF file
    using the Scientific.IO.NetCDF.NetCDFFile interface.
    """
    __metaclass__ = classreg
    @classmethod
    def isMine(cls, *args, **kwds):
        try:
            cls(*args, **kwds)
            return True
        except:
            return False
            
    def __new__(cls, *args, **kwds):
        new = object.__new__(cls, *args, **kwds)
        new.variables = OrderedDict()
        new.dimensions = OrderedDict()
        new._ncattrs = ()
        new._operator_exclude_vars = ()
        return new
    
    def __init__(self, *args, **properties):
        for k, v in properties.iteritems():
            setattr(self, k, v)

    def __setattr__(self, k, v):
        if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
            if k not in self._ncattrs:
                self._ncattrs += (k, )
        object.__setattr__(self, k, v)
    def __delattr__(self, k):
        if k in self._ncattrs:
            self._ncattrs = tuple([k_ for k_ in self._ncattrs if k_ != k])
        object.__delattr__(self, k)

    def createDimension(self, name, length):
        """
        name - string name for dimension
        length - maximum length of dimension
        """
        dim = self.dimensions[name] = PseudoNetCDFDimension(self, name, length)
        return dim

    def createVariable(self, name, type, dimensions, **properties):
        """
        name - string
        type - numpy dtype code (e.g., 'f', 'i', 'd')
        dimensions - tuple of dimension keys that can be
                     found in objects' dimensions dictionary
        """
        import numpy as np
        if type == 'S': type = 'c'
        if isinstance(properties.get('values', 1), np.ma.MaskedArray) or 'fill_value' in properties:
            var = self.variables[name] = PseudoNetCDFMaskedVariable(self, name, type, dimensions, **properties)
        else:
            var = self.variables[name] = PseudoNetCDFVariable(self, name, type, dimensions, **properties)
        return var

    def close(self):
        """
        Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
        """
        pass

    def ncattrs(self):
        return self._ncattrs

    def __add__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '+', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __sub__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '-', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __mul__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '*', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __div__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '/', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __floordiv__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '//', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __pow__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '**', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __and__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '&', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __or__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '|', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)

    def __xor__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '^', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)

    def __mod__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '%', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __lt__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '<', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __gt__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '>', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)

    def __eq__(self, lhs):
        if isinstance(lhs, (NetCDFFile, PseudoNetCDFFile)):
            from _functions import pncbo
            return pncbo(op = ' == ', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
        else:
            return lhs.__eq__(self)

    def __le__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '<=', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    def __ge__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '>=', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)

    def __ne__(self, lhs):
        from _functions import pncbo
        return pncbo(op = '!=', ifile1 = self, ifile2 = lhs, verbose = False, coordkeys = self._operator_exclude_vars)
    
    
    sync = close
    flush = close

class PseudoNetCDFFileMemmap(PseudoNetCDFFile):
    """
    Provides basic PseudoNetCDFFile functionality, but
    does not require that variables be created in memmory
    """
    def createVariable(self, name, type, dimensions, map, keep = True):
        var = PseudoNetCDFVariableMemmap(self, name, type, dimensions, map)
        if keep:
            self.variables[name] = var
        return var

class PseudoNetCDFVariables(OrderedDefaultDict):
    """
    PseudoNetCDFVariables provides a special implementation
    of the default dictionary that provides efficient access 
    to variables of a PseudoNetCDFFile.  PseudoNetCDFFiles may
    have large variables that should only be loaded if accessed.
    PseudoNetCDFVariables allows a user to specify a function
    that can create variables on demand.
    """
    def __init__(self, func, keys):
        """
        func: Function that takes a key and provides a 
              PseudoNetCDFVariable
        keys: list of keys that the dictionary should
              act as if it has
        """
        super(PseudoNetCDFVariables, self).__init__()
        self.__func = func
        self.__keys = keys
    def __missing__(self, k):
        """
        If the dictionary does not have a key, check if the
        user has provided that key.  If so, call the user 
        specifie function to create the variable.
        """
        if k in self.keys():
            return self.__func(k)
        else:
            raise KeyError, 'missing "%s"' % (k, )

    def addkey(self, k):
        """
        Allow the user to extend keys after the object
        has been created.
        """
        self.__keys.append(k)

    def keys(self):
        return tuple(set(dict.keys(self) + self.__keys))
    
    def has_key(self, k):
        return k in self.keys()
    
    def iteritems(self):
        for k in self.keys():
            yield k, self[k]
    
    def iterkeys(self):
        for k in self.keys():
            yield k

class PseudoNetCDFTest(unittest.TestCase):
    def setUp(self):
        self.tncf = PseudoNetCDFFile()
        
    def testNetCDFFileInit(self):
        from numpy import arange
        tncf = self.tncf
        self.assert_(tncf.variables == {})
        self.assert_(tncf.dimensions == {})

        tncf.createDimension('TIME', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        self.assert_(len(tncf.dimensions['TIME']) == 24)
        self.assert_(len(tncf.dimensions['LAY']) == 4)
        self.assert_(len(tncf.dimensions['ROW']) == 5)
        self.assert_(len(tncf.dimensions['COL']) == 6)
        
        tncf.fish = 2
        setattr(tncf, 'FROG-DOG', 'HAPPY')

        o3 = tncf.createVariable('O3', 'f', ('TIME', 'LAY', 'ROW', 'COL'))
        self.assert_(tncf.variables.keys() == ['O3'])
        
        o3[:] = arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)
        o3.units = 'ppbv'
        self.assert_((o3 == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        self.assert_((tncf.variables['O3'] == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        self.assert_(o3.typecode() == 'f')

        filedims = tncf.dimensions.keys()
        filedims.sort()
        vardims = list(o3.dimensions)
        vardims.sort()

        self.assert_(filedims == vardims)
        from PseudoNetCDF.pncgen import Pseudo2NetCDF
        n = Pseudo2NetCDF().convert(tncf)
        self.assertEqual(n.variables.keys(), ['O3'])
        self.assertEqual(dict([(k, len(v)) for k, v in n.dimensions.iteritems()]), {'TIME': 24, 'LAY': 4, 'ROW': 5, 'COL': 6})
        self.assert_((n.variables['O3'][...] == tncf.variables['O3'][...]).all())
        self.assert_(n.variables['O3'].units == 'ppbv')
        self.assertEqual(n.fish, 2)
        self.assertEqual(getattr(n, 'FROG-DOG'), 'HAPPY')
        
    def testNetCDFVariables(self):
        from numpy import arange
        tncf = PseudoNetCDFFile()
        tncf.createDimension('TIME', 24)
        tncf.createDimension('LAY', 4)
        tncf.createDimension('ROW', 5)
        tncf.createDimension('COL', 6)

        const = lambda *args, **kwds: PseudoNetCDFVariable(tncf, args[0], 'f', ('TIME', 'LAY', 'ROW', 'COL'), values = arange(24 * 4 * 5 * 6).reshape((24, 4, 5, 6)))
        tncf.variables = PseudoNetCDFVariables(const, ['NO', 'O3'])
        self.assert_((tncf.variables['O3'] == arange(24 * 4 * 5 * 6).reshape(24, 4, 5, 6)).all())
        
        
    def runTest(self):
        pass

if __name__ == '__main__':
    unittest.main()
