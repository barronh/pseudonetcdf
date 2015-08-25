__doc__ = r"""
.. _dumper
:mod:`dumper` -- PseudoNetCDF dump module
============================================

.. module:: pncdump
   :platform: Unix, Windows
   :synopsis: Provides ncdump like functionaility for PseudoNetCDF
.. moduleauthor:: Barron Henderson <barronh@unc.edu>

"""

__all__=['pncload',]

import code
import atexit
import os

class PNCConsole(code.InteractiveConsole):
    """
    The PNCConsole is designed to create a fast Integrated
    Development Environment for scientific analysis.
    """
    def __init__(self, locals = None, filename = '<console>',
                histfile = os.path.expanduser("~/.pncload-history")):
        """
        locals - dictionary of local variables
        filename - filename for error printouts
        histfile - path for history
        """
        code.InteractiveConsole.__init__(self, locals = locals, filename = filename)
        self.init_history(histfile)
    
    def init_history(self, histfile):
        """
        Prepare history for use from histfile (typically last session)
        """
        import readline
        readline.parse_and_bind("tab: complete")
        if hasattr(readline, "read_history_file"):
            try:
                readline.read_history_file(histfile)
                readline.set_history_length(1000)
            except IOError:
                pass
            atexit.register(self.save_history, histfile)
            
    def save_history(self, histfile):
        """
        Write history from session to disk
        """
        import readline
        readline.write_history_file(histfile)

_torem = ['(', ')', '[', ']', '{', '}', '@', ',', ':', '.', '`', '=', ';', '+=', '-=', '*=', '/=', '//=', '%=', '&=', '|=', '^=', '>>=', '<<=', '**=', '+', '-', '*', '**', '/', '//', '%', '<<', '>>', '&', '|', '^', '~', '<', '>', '<=', '>=', '==', '!=', '<>']
def _clean(p):
    """
    Remove operator symbols from file names
    """
    for d in _torem:
        p = os.path.basename(p).replace(d, '')
    return p

def createconsole(ifiles, options):
    """
    Use standard pncparse ifiles and options to create
    a working environment.
    
    1) Access files by a short name or indexed name (ifile%d)
    2) Access variables by name and file index (var_%d)
    3) Has access to pylab and pncview functions
    """
    console = PNCConsole()
    exec("from pylab import *", None, console.locals)
    exec("from PseudoNetCDF.pncview import *; interactive(True)", None, console.locals)
    ipaths = [_clean(ipath) for ipath in options.ipath]
    spaths = [ipath[:6] for ipath in ipaths]
    spathsoc = dict([(k, 0) for k in spaths])
    spathso = []
    npathso = []
    for filei, spath in enumerate(spaths):
        spathso.append(spath + '_' + str(spathsoc[spath]))
        spathsoc[spath] += 1
        npathso.append('ifile%d' % filei)
    for rpath, npath, spath in zip(options.ipath, npathso, spathso):
        print spath + ' = ' + npath + ' = ' + rpath
    for filei, (ipath, npath, spath, ifile) in enumerate(zip(ipaths, npathso, spathso, ifiles)):
        console.locals[npath] = ifile
        try:
            exec(ipath + ' = ' + npath, None, console.locals)
        except: pass
        try:
            exec(spath + ' = ' + npath, None, console.locals)
        except: pass
        console.locals.update(dict([('%s_%d' % (k, filei), v) for k, v in ifile.variables.iteritems()]))
    return console

def main():
    from pncparse import pncparse
    ifiles, options = pncparse(has_ofile = False, interactive = False)
    console = createconsole(ifiles, options)
    console.interact()
    
if __name__ == '__main__':
    main()