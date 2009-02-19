__all__ = ['cartesian', 'sliceit']
__doc__ = """
.. _util
:mod:`util` -- CAMx basic util
==============================

.. module:: util
   :platform: Unix, Windows
   :synopsis: Provides simple utilites for camxfiles
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

def cartesian(x, y):
    """Iterator for an 'outer' or cartesian join of
    x iterator and y iterator
    """
    for i in x:
        for j in y:
            yield i,j

def sliceit(args):
    """
    If arguments are a slice, return it.  Otherwise,
    create a slice.
    """
    if isinstance(args,slice):
       return args
       
    try:
        return slice(*args)
    except TypeError:
        return slice(args,args+1)

