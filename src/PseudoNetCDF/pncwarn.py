import sys
import warnings

std_showwarning = warnings.showwarning


def warn(*args, **kwds):
    """
    Thin wrapper around warnings.warn to prepend **PNC:
    to warnings. This makes it easy to find PseudoNetCDF
    warnings.

    Arguments
    ---------
    args : arguments for warnings.warn
    kwds : keywords for warnings.warn

    Returns
    -------
    out : same as warnings.warn
    """

    warnings.showwarning = clean_showwarning
    out = warnings.warn(*args, **kwds)
    warnings.showwarning = std_showwarning
    return out


def clean_showwarning(message, category, filename, lineno, file=None,
                      line=None):
    global _first_read_only
    if file is None:
        file = sys.stderr
        if file is None:
            # sys.stderr is None when run with pythonw.exe - warnings get lost
            return
    try:
        warnstr = '**PNC:%s:%s:%s:\n  %s\n' % (
            (filename), lineno, category.__name__, message)
        file.write(warnstr)
    except OSError:
        pass  # the file (probably stderr) is invalid - this warning gets lost.
    return


filterwarnings = warnings.filterwarnings
simplefilter = warnings.simplefilter


def quiet():
    simplefilter('ignore', UserWarning)
