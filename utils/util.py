__doc__ = """
.. _util
:mod:`util` -- pyPA utilities
==============================

.. module:: util
   :platform: Unix, Windows
   :synopsis: Provides unclassified general utilities used in :ref:`pyPA`
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__ = ['defaultdict', 'getheights', 'stopwatch', 'descendAttrDict', 'AttrDict', 'add_tuple', 'toints', 'superTuple', 'AlwaysEquals']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from warnings import warn
from datetime import datetime

class defaultdict(dict):
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and                 not hasattr(default_factory, '__call__')):
            raise TypeError('first argument must be callable')
        dict.__init__(self, *a, **kw)
        self.default_factory = default_factory
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value
    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()
    def copy(self):
        return self.__copy__()
    def __copy__(self):
        return type(self)(self.default_factory, self)
    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))
    def __repr__(self):
        return 'defaultdict(%s, %s)' % (self.default_factory,
                                        dict.__repr__(self))        

def getheights(file):
    warn("get heights has been deprecated in favor of hp.variables['HGHT']",DeprecationWarning)
    from pyPA.utils.CAMxFiles import height_pressure
    import operator
    hp=height_pressure(file)
    heights=[]
    for d,t,k in hp:
        thisht=hp.seekandread(d,t,k)
        heights.append(reduce(operator.add,thisht)/len(thisht))
    return heights

class stopwatch:
    """
    Example use
    t=stopwatch()
    print 1
    print t()
    """
    def __init__(self):
        self.start=datetime.now()
    def __call__(self):
        new=datetime.now()
        result=new-self.start
        self.start=new
        return result

class AttrDict(dict):
    """For easy access to values"""
    def __new__(cls,*args,**kwds):
        result = dict.__new__(cls,*args,**kwds)
        return result
        
    def __init__(self,*args,**kwds):
        dict.__init__(self, *args, **kwds)
        for k,v in self.iteritems():
            if type(v)==dict:
                self[k]=AttrDict(v)
            else:
                self[k]=v
                
    def __setattr__(self, attr, value):
        self[attr] = value
        
    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError(attr)

descendAttrDict=AttrDict


def add_tuple(t1, t2):
    """Add tuple 1 values to tuple 2 values
    """
    return tuple([i + j
                  for i, j in zip(t1, t2)])

def toints(line,delim=None):
    """Convert tuple values to ints
    """
    return tuple([int(i)
                  for i in line.split(delim)])

def superTuple(name, attributes):
    """Creates a Super Tuple class.
    From http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/218485
    """
    dct = {}
    #Create __new__.
    nargs = len(attributes)
    def _new_(cls, *args):
        if len(args) != nargs:
            raise TypeError("%s takes %d arguments (%d given)." % (cls.__name__,
                                                                   nargs,
                                                                   len(args)))
        return tuple.__new__(cls, args)
    dct["__new__"] = staticmethod(_new_)
    #Create __repr__.
    def _repr_(self):
        contents = [repr(elem) for elem in self]
        return "%s<%s>" % (self.__class__.__name__,
                           ", ".join(contents))
    dct["__repr__"] = _repr_
    #Create attribute properties.
    def getter(i):
        return lambda self: self.__getitem__(i)
    for index, attribute in enumerate(attributes):
        dct[attribute] = property(getter(index))
    #Set slots.
    dct["__slots__"] = []
    #Return class.
    return type(name, (tuple,), dct)


class AlwaysEquals:
    """Will always return true...
    """
    def __eq__(self, other):
        return True
    
