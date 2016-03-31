#!/usr/bin/env python
# Run this script with pnc options
from __future__ import print_function

import sys
import os

from PseudoNetCDF.pncload import PNCConsole
from PseudoNetCDF.pncparse import pncparse
from PseudoNetCDF.coordutil import gettimes
import warnings
warn=warnings.warn
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
Normalize = matplotlib.colors.Normalize

LogNorm = matplotlib.colors.LogNorm
SymLogNorm = matplotlib.colors.SymLogNorm
BoundaryNorm = matplotlib.colors.BoundaryNorm

LogFormatter = matplotlib.ticker.LogFormatter
ScalarFormatter = matplotlib.ticker.ScalarFormatter

fig = plt.figure()
ax = fig.add_subplot(111)
def make1d(ifile, options):
    for fi, ifile in enumerate(ifiles):
        variables = options.variables
        if variables is None:
            variables = [key for key, var in ifile.variables.items() if var.ndim == 2]
        if len(variables) == 0:
            raise ValueError('Unable to heuristically determin plottable variables; use -v to specify variables for plotting')
        for varkey in variables:
            if not options.overlay:
                ax.cla()
            var = ifile.variables[varkey]
            vals = var[:]
            label = getattr(var, 'standard_name', varkey).strip()
            varunit = getattr(var, 'units', 'unknown').strip()
            print(varkey, sep = '')
            dimkey = var.dimensions[0]
            ax.set_xlabel(dimkey)
            if options.time:
                x = gettimes(ifile)
            elif dimkey in ifile.variables:
                x = ifile.variables[var.dimensions[0]][:]
            else:
                x = np.arange(var.shape[0])
            patches = ax.plot(x, vals, label = label)
            ax.set_ylabel(varunit) 
            plt.legend()
            fmt = 'png'
            outpath = options.outpath
            if len(ifiles) > 1:                
                outpath += ('_%%0%dd' % len(str(len(ifiles)))) % fi
            figpath = os.path.join(outpath + varkey + '.' + fmt)
            if options.interactive:
                csl = PNCConsole(locals = globals())
                csl.interact()
            
            fig.savefig(figpath)
            print('Saved fig', figpath)
        
if __name__ == '__main__':
    from PseudoNetCDF.pncparse import pncparse, getparser
    parser = getparser(has_ofile = True, plot_options = True, interactive = True)
    parser.add_argument('--time', action = 'store_true', help = 'Use time object as x-axis')
    ifiles, options = pncparse(has_ofile = True, plot_options = True, interactive = True, parser = parser)
    make1d(ifiles, options)
