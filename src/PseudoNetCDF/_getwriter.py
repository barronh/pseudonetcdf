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
    """See PseudoNetCDF.pncgen.pncgen

    *args : iterable
    **kwds : keywords
        keywords for pncgen
    help : boolean
        without format, returns help of pncopen and with format keyword,
        returns help of the function that writes that format.
    """
    from PseudoNetCDF.pncgen import pncgen
    formathelp = kwds.pop('help', False)
    if formathelp:
        if 'format' in kwds:
            format = kwds['format']
            writer = getwriterdict()[format]
            return help(writer)
        else:
            return help(pncgen)
    return pncgen(*args, **kwds)
