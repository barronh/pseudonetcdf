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

def main():
    from pncparse import pncparser
    ifiles, options = pncparser(has_ofile = False)
    from permm.Shell import PERMConsole
    console = PNCConsole()
    exec("from pylab import *", None, console.locals)
    exec("from PseudoNetCDF.pncview import *; interactive(True)", None, console.locals)
    for filei, ifile in enumerate(ifiles):
        console.locals['ifile%d' % filei] = ifile
    console.interact()
    
if __name__ == '__main__':
    main()