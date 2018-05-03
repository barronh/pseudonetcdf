__all__ = ['getwriterdict', 'registerwriter']

_writers = []


def testwriter(writer, *args, **kwds):
    try:
        writer(*args, **kwds)
        return True
    except Exception:
        return False


def registerwriter(name, writer):
    global _writers
    _writers.insert(0, (name, writer))


def getwriterdict():
    return dict(_writers)


def pncwrite(*args, **kwds):
    """See PseudoNetCDF.pncgen.pncgen"""
    from PseudoNetCDF.pncgen import pncgen
    return pncgen(*args, **kwds)
