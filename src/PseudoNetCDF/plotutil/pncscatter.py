__all__ = ['pncscatter']

import matplotlib.pyplot as plt
import numpy as np


def pncscatter(args):
    ifiles = args.ifiles
    fig = plt.figure()
    if len(args.figure_keywords) > 0:
        plt.setp(fig, **args.figure_keywords)
    sax = fig.add_subplot(111)  # add_axes([.1, .15, .8, .8])
    if len(args.axes_keywords) > 0:
        plt.setp(sax, **args.axes_keywords)

    ifilex = ifiles[0]
    ifiley = ifiles[1]

    for target in args.variables:
        sax.set_xlabel('Time (UTC)')
        varx = ifilex.variables[target]
        varxdesc = getattr(varx, 'description', None)
        vary = ifiley.variables[target]
        varydesc = getattr(vary, 'description', None)
        sax.set_ylabel(varydesc)
        sax.set_xlabel(varxdesc)
        del sax.lines[:]
        vmin = np.minimum(varx[:].min(), vary[:].min())
        vmax = np.maximum(varx[:].max(), vary[:].max())
        sax.plot([vmin, vmax], [vmin, vmax], color='k')
        sax.plot(varx[:].ravel(), vary[:].ravel(),
                 marker='o', ls='none', markeredgecolor='none')
        sax.set_xlim(vmin, vmax)
        sax.set_ylim(vmin, vmax)
        # plt.setp(sax.xaxis.get_ticklabels(),rotation = 45)
        figpath = args.outpath + target + '.' + args.figformat
        for pc in args.plotcommands:
            exec(pc)
        fig.savefig(figpath)
        if args.verbose > 0:
            print('Saved fig', figpath)
