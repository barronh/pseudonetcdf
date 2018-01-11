import sys
import warnings
warn = warnings.warn
def clean_showwarning(message, category, filename, lineno, file = None, line = None):
    global _first_read_only
    if file is None:
        file = sys.stderr
        if file is None:
            # sys.stderr is None when run with pythonw.exe - warnings get lost
            return
    try:
        warnstr = '**PNC:%s:%s:%s:\n  %s\n' % ((filename), lineno, category.__name__, message)
        file.write(warnstr)
    except OSError:
        pass # the file (probably stderr) is invalid - this warning gets lost.
    return
warnings.showwarning = clean_showwarning

filterwarnings = warnings.filterwarnings
simplefilter = warnings.simplefilter

def quiet() :
    simplefilter('ignore', UserWarning)

