#!/usr/bin/env python
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from PseudoNetCDF.coordutil import gettimes
import numpy as np
import os

def plotts(args):
    ifiles = args.ifiles
    times = [gettimes(ifile) for ifile in ifiles]
    fig = plt.figure()
    ax = fig.add_subplot(111)#add_axes([.1, .15, .8, .8])
    ax.set_xlabel('Time (UTC)')
    split = 25
    for target in args.variables:
        vars = [ifile.variables[target] for ifile in ifiles]
        unit = getattr(vars[0], 'units', 'unknown')
        ax.set_ylabel(target + '(' + unit + ')')
        del ax.lines[:]
        nvars = len(vars)
        varwidth = .8/nvars/1.1
        po = np.arange(24) + 0.1 + varwidth/2
        for vi, (time, var) in enumerate(zip(times, vars)):
            vals = var[:]
            if args.squeeze:
                vals = vals.squeeze()
            vardesc = getattr(var, 'description', None)
            varb = ax.plot(time, vals[:], label = vardesc)
        #plt.setp(ax.xaxis.get_ticklabels(),rotation = 45)
        plt.legend()
        figpath = args.outpath + target + '.' + args.figformat
        for pc in args.plotcommands:
           exec(pc, globals(), locals())
        fig.savefig(figpath)
        if args.verbose > 0: print('Saved fig', figpath)
        print(figpath)

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
    plotts(args)
