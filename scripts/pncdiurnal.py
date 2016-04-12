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
    sax.set_xlabel('Time (UTC)')
    split = 25
    for target in args.variables:
        vars = [ifile.variables[target] for ifile in ifiles]
        unit = getattr(vars[0], 'units', 'unknown')
        sax.set_ylabel(target + '(' + unit + ')')
        vals = [var[:] for var in vars]
        hvars = [[np.ma.compressed(val[hour == i]) for i in range(24)] for hour, val in zip(hours, vals)]
        del sax.lines[:]
        nvars = len(vars)
        varwidth = .8/nvars/1.1
        po = np.arange(24) + 0.1 + varwidth/2
        for vi, var in enumerate(hvars):
            varb = sax.boxplot(var, positions = po + vi * varwidth * 1.1, widths = varwidth, patch_artist = True)
            color = plt.get_cmap()((vi + .5)/float(nvars))
            plt.setp(varb.values(), color = color)
            plt.setp(varb['medians'], color = 'k')
            plt.setp(varb['fliers'], markeredgecolor = color)
        sax.set_xlim(-.5, 24.5)
        sax.set_xticks(range(0, 25))
        sax.set_xticklabels([str(i) for i in range(0, 25)])
        #plt.setp(sax.xaxis.get_ticklabels(),rotation = 45)
        figpath = args.outpath + target + '_box.png'
        fig.savefig(figpath)
        print figpath

if __name__ == '__main__':
    from PseudoNetCDF.pncparse import getparser, pncparse
    parser = getparser(plot_options = True, has_ofile = True)
    parser.epilog += """
    -----
box.py inobs inmod target [target ...]
inobs - path to obs
inmod - path to mod
target - variable name
"""
    ifiles, args = pncparse(plot_options = True, has_ofile = True, parser = parser)
    plot_diurnal_box(ifiles, args)
