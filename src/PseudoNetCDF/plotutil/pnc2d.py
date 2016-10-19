#!/usr/bin/env python
# Run this script with pnc options
from __future__ import print_function

import sys
import os

from PseudoNetCDF.pncload import PNCConsole
from PseudoNetCDF.pncparse import pncparse
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

def make2ds(args):
    ifiles = args.ifiles
    if len(args.figure_keywords) > 0:
        plt.setp(fig, **args.figure_keywords)
    if len(args.axes_keywords) > 0:
        plt.setp(ax, **args.axes_keywords)
    nborders = len(ax.collections)
    for fi, ifile in enumerate(ifiles):
        variables = args.variables
        if variables is None:
            variables = [key for key, var in ifile.variables.items() if var.ndim == 2]
        if len(variables) == 0:
            raise ValueError('Unable to heuristically determin plottable variables; use -v to specify variables for plotting')
        for varkey in variables:
            var = ifile.variables[varkey]
            vals = var[:]
            if args.squeeze:
                vals = vals.squeeze()
            
            if args.normalize is None:
                from scipy.stats import normaltest
                vmin, vmax = vals.min(), vals.max()
                if normaltest(vals.ravel())[1] < 0.05:
                    cvals = np.ma.compressed(vals)
                    boundaries = np.percentile(cvals, np.arange(0, 110, 10))
                    warn('Autoselect deciles colormap of %s; override width --norm' % varkey)
                else:
                    boundaries = np.linspace(vmin, vmax, num = 11)
                    warn('Autoselect linear colormap of %s; override width --norm' % varkey)
                if (boundaries.max() / np.ma.masked_values(boundaries, 0).min()) > 10000:
                    formatter = LogFormatter(labelOnlyBase = False)
                else:
                    formatter = None
                norm = BoundaryNorm(boundaries, ncolors = 256)
            else:
                norm = eval(args.normalize)
                formatter = None
            if not args.colorbarformatter is None:
                try:
                    formatter = eval(args.colorbarformatter)
                except:
                    formatter = args.colorbarformatter

                
            vmin, vmax = vals.min(), vals.max()
            if not norm.vmin is None:
                vmin = norm.vmin
            if not norm.vmax is None:
                vmax = norm.vmax
            
            varunit = getattr(var, 'units', 'unknown').strip()
            vardims = [dk for dk, dv in zip(var.dimensions, var.shape) if dv != 1]
            print(varkey, sep = '')
            del ax.collections[nborders:]
            if args.swapaxes:
                patches = ax.pcolor(vals.T, norm = norm)
                ax.set_xlabel(vardims[0])
                ax.set_ylabel(vardims[1]) 
            else:
                patches = ax.pcolor(vals, norm = norm)
                ax.set_xlabel(vardims[1])
                ax.set_ylabel(vardims[0])

            height = vals.shape[0]
            width = vals.shape[1]
            if width >= height:
                orientation = 'horizontal'
            else:
                orientation = 'vertical'
            try:
                cax = cbar.ax
                cax.cla()
            except:
                cax = None
            if vals.max() > vmax and vals.min() < vmin:
                extend = 'both'
            elif vals.max() > vmax:
                extend = 'max'
            elif vals.min() < vmin:
                extend = 'min'
            else:
                extend = 'neither'
            cbar = fig.colorbar(patches, orientation = orientation, cax = cax, extend = extend, format = formatter)
            del cbar.ax.texts[:]
            cbar.set_label(varkey + ' (' + varunit + '; min=%.3g; max=%.3g)' % (var[:].min(), var[:].max()))
 #           if orientation == 'vertical':
 #               cbar.ax.text(.5, 1.05, '%.3g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
#                cbar.ax.text(.5, -.06, '%.3g ' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
#            else:
#                cbar.ax.text(1.05, .5, ' %.3g' % var[:].max(), verticalalignment = 'center', horizontalalignment = 'left')
#                cbar.ax.text(-.06, .5, '%.3g ' % var[:].min(), verticalalignment = 'center', horizontalalignment = 'right')
            #cbar.update_ticks()
            fmt = 'png'
            outpath = args.outpath
            if len(ifiles) > 1:                
                lstr = str(fi).rjust(len(str(len(ifiles))), '0')
            else:
                lstr = ''
                
            figpath = os.path.join(outpath + varkey + lstr + '.' + fmt)
            if args.interactive:
                csl = PNCConsole(locals = globals())
                csl.interact()
            
            fig.savefig(figpath)
            if args.verbose > 0: print('Saved fig', figpath)
        
if __name__ == '__main__':
    from PseudoNetCDF.pncparse import pncparse, getparser
    parser = getparser(has_ofile = True, plot_options = True, interactive = True)
    #parser.add_argument('--no-squeeze', dest = 'squeeze', default = True, action = 'store_false', help = 'Squeeze automatically removes singleton dimensions; disabling requires user to remove singleton dimensions with --remove-singleton option')
    parser.add_argument('--swapaxes', action = 'store_true', help = 'Swap x-y axes')
    ifiles, args = pncparse(has_ofile = True, plot_options = True, interactive = True, parser = parser)
    make2d(ifiles, args)
