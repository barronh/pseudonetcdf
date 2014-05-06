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
import readline
import atexit
import os

class PNCConsole(code.InteractiveConsole):
    def __init__(self, locals = None, filename = '<console>',
                histfile = os.path.expanduser("~/.pncload-history")):
        code.InteractiveConsole.__init__(self)
        self.init_history(histfile)
    
    def init_history(self, histfile):
        readline.parse_and_bind("tab: complete")
        if hasattr(readline, "read_history_file"):
            try:
                readline.read_history_file(histfile)
                readline.set_history_length(1000)
            except IOError:
                pass
            atexit.register(self.save_history, histfile)
            
    def save_history(self, histfile):
        readline.write_history_file(histfile)

_torem = ['(', ')', '[', ']', '{', '}', '@', ',', ':', '.', '`', '=', ';', '+=', '-=', '*=', '/=', '//=', '%=', '&=', '|=', '^=', '>>=', '<<=', '**=', '+', '-', '*', '**', '/', '//', '%', '<<', '>>', '&', '|', '^', '~', '<', '>', '<=', '>=', '==', '!=', '<>']
def _clean(p):
    for d in _torem:
        p = p.replace(d, '')
    return p

def createconsole(ifiles, options):
    console = PNCConsole()
    exec("from pylab import *", None, console.locals)
    exec("from PseudoNetCDF.pncview import *; interactive(True)", None, console.locals)
    ipaths = [_clean(ipath) for ipath in options.ifile]
    spaths = [ipath[:6] for ipath in ipaths]
    spathsoc = dict([(k, 0) for k in spaths])
    spathso = []
    for spath in spaths:
        spathso.append(spath + '_' + str(spathsoc[spath]))
        spathsoc[spath] += 1
    
    for filei, (ipath, spath, ifile) in enumerate(zip(ipaths, spathso, ifiles)):
        npath = 'ifile%d' % filei
        console.locals[npath] = ifile
        exec(ipath + ' = ' + npath, None, console.locals)
        exec(spath + ' = ' + npath, None, console.locals)
        console.locals.update(dict([('%s_%d' % (k, filei), v) for k, v in ifile.variables.iteritems()]))
    return console

def main():
    from pncparse import pncparser
    ifiles, options = pncparser(has_ofile = False)
    console = createconsole(ifiles, options)
    console.interact()
    
if __name__ == '__main__':
    main()