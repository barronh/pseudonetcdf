#!/usr/bin/env python
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from PseudoNetCDF.coordutil import gettimes
import numpy as np
import os

def plot_diurnal_box(ifiles, args):
    times = [gettimes(ifile) for ifile in ifiles]
    hours = [np.array([t.hour for t in time]) for time in times]
    fig = plt.figure()
    sax = fig.add_subplot(111)
    if len(args.figure_keywords) > 0:
        plt.setp(fig, **args.figure_keywords)
    if len(args.axes_keywords) > 0:
        plt.setp(sax, **args.axes_keywords)
    sax.set_xlabel('Time (UTC)')
    split = 25
    for target in args.variables:
        vars = [ifile.variables[target] for ifile in ifiles]
        unit = getattr(vars[0], 'units', 'unknown')
        sax.set_ylabel(target + ' (' + unit + ')')
        vals = [var[:] for var in vars]
        hvars = [[np.ma.compressed(val[hour == i]) for i in range(24)] for hour, val in zip(hours, vals)]
        del sax.lines[:]
        nvars = len(vars)
        varwidth = .8/nvars/1.1
        po = np.arange(24) + 0.1 + varwidth/2
        ncolors = max(float(nvars), 10)
        try:
            from cycler import cycler
            props = cycler('color', [plt.get_cmap()((nc + .5)/float(ncolors)) for nc in range(ncolors)])
        except:
            props = [dict(color = plt.get_cmap()((nc + .5)/float(ncolors))) for nc in range(ncolors)]
        for vi, (var, propd) in enumerate(zip(hvars, props)):
            varb = sax.boxplot(var, positions = po + vi * varwidth * 1.1, widths = varwidth, patch_artist = True)
            plt.setp([i for i in varb.values()], **propd)
            plt.setp(varb['medians'], color = 'k')
            plt.setp(varb['fliers'], markeredgecolor = propd['color'])
        sax.set_xlim(-.5, 24.5)
        sax.set_xticks(range(0, 25))
        sax.set_xticklabels([str(i) for i in range(0, 25)])
        #plt.setp(sax.xaxis.get_ticklabels(),rotation = 45)
        figpath = args.outpath + target + '.' + args.figformat
        fig.savefig(figpath)
        if args.verbose > 0: print('Saved fig', figpath)

if __name__ == '__main__':
    from PseudoNetCDF.pncparse import getparser, pncparse
    parser = getparser(plot_options = True, has_ofile = True)
    parser.add_argument('--whis', default = None, help = 'See pydoc matplotlib.axes.Axes.boxplot')
    parser.epilog += """
    -----
box.py inobs inmod target [target ...]
inobs - path to obs
inmod - path to mod
target - variable name
"""
    ifiles, args = pncparse(plot_options = True, has_ofile = True, parser = parser)
    plot_diurnal_box(ifiles, args)
