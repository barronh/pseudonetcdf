#!/usr/bin/env python
import sys
import numpy as np
from warnings import warn
from collections import defaultdict
from netCDF4 import MFDataset, Dataset
from datetime import datetime, timedelta
from PseudoNetCDF.coordutil import getmap


def plot(ifiles, args):
    import pylab as pl
    from pylab import figure, NullFormatter, close, rcParams
    from PseudoNetCDF.coordutil import getsigmamid, getpresmid, getpresbnds, getsigmabnds
    rcParams['text.usetex'] = False
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, LogNorm
    map = not args.nomap
    scale = args.scale
    minmax = eval(args.minmax)
    minmaxq = eval(args.minmaxq)
    sigma = args.sigma
    maskzeros = args.maskzeros
    try:
        f, = ifiles
    except:
        raise ValueError('curtain plot expects one file when done. Try stack time --stack=time to concatenate')
    
    if sigma:
        vertcrd = getsigmabnds(f)
    else:
        vertcrd = getpresbnds(f, pref = 101325., ptop = getattr(f, 'VGTOP', 10000))
        if vertcrd.max() > 2000:  vertcrd /= 100.


    reversevert = not (np.diff(vertcrd) > 0).all()
    for var_name in args.variables:
        temp = defaultdict(lambda: 1)
        try:
            eval(var_name, None, temp)
            var = eval(var_name, None, f.variables)[:]
        except:
            temp[var_name]
            var = f.variables[var_name][:]
        if args.itertime:
            vars = [('time%02d' % vi, v) for vi, v in enumerate(var)]
        else:
            if var.shape[0] != 1:
                parser.print_help()
                sys.stderr.write('\n*****\n\nFile must have only one time or use the itertime options; to reduce file to one time, see the --slice and --reduce operators\n\n')
                exit()
            vars = [('', var[0])]
        
        for lstr, var in vars:
            bmap = None
            if maskzeros: var = np.ma.masked_values(var, 0)
            vmin, vmax = np.percentile(np.ma.compressed(var).ravel(), list(minmaxq))
            if minmax[0] is not None:
                vmin = minmax[0]
            if minmax[1] is not None:
                vmax = minmax[1]
            
            if not args.normalize is None:
                norm = eval(args.normalize)
                if norm.scaled():
                    vmin = norm.vmin; vmax = norm.vmax
            else:
                if scale == 'log':
                    bins = np.logspace(np.log10(vmin), np.log10(vmax), 11)
                elif scale == 'linear':
                    bins = np.linspace(vmin, vmax, 11)
                elif scale == 'deciles':
                    bins = np.percentile(np.ma.compressed(np.ma.masked_greater(np.ma.masked_less(var, vmin), vmax)).ravel(), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
                    bins[0] = vmin; bins[-1] = vmax
                norm = BoundaryNorm(bins, ncolors = 256)
            
            if map:
                fig = pl.figure(figsize = (8, 8))
                fig.subplots_adjust(hspace = .3, wspace = .3)
                axmap = fig.add_subplot(3,3,5)
                try:
                    cmaqmap = getmap(f)
                    cmaqmap.drawcoastlines(ax = axmap)
                    cmaqmap.drawcountries(ax = axmap)
                    if args.states: cmaqmap.drawstates(ax = axmap)
                    if args.counties: cmaqmap.drawcounties(ax = axmap)
                    cmaqmap.drawparallels(np.arange(-90, 100, 10), labels = [True, True, False, False], ax = axmap)
                    cmaqmap.drawmeridians(np.arange(-180, 190, 20), labels = [False, False, True, True], ax = axmap)
                except Exception as e:
                    warn('An error occurred and no map will be shown:\n%s' % str(e))
                axn = fig.add_subplot(3,3,2, sharex = axmap)
                axw = fig.add_subplot(3,3,4, sharey = axmap)
                axe = fig.add_subplot(3,3,6, sharey = axmap)
                axs = fig.add_subplot(3,3,8, sharex = axmap)
                cax = fig.add_axes([.8, .7, .05, .25])
                for ax in [axmap, axe]:
                    ax.yaxis.set_major_formatter(NullFormatter())
                for ax in [axmap, axn]:
                    ax.xaxis.set_major_formatter(NullFormatter())
                for ax in [axn, axs]:
                    if sigma:
                        ax.set_ylabel('sigma')
                    else:
                        ax.set_ylabel('pressure')
                for ax in [axe, axw]:
                    if sigma:
                        ax.set_xlabel('sigma')
                    else:
                        ax.set_xlabel('pressure')
                xyfactor = 1
            else:
                fig = pl.figure(figsize = (16, 4))
                fig.subplots_adjust(bottom=0.15)
                axw = fig.add_subplot(1,4,1)
                axn = fig.add_subplot(1,4,2)
                axe = fig.add_subplot(1,4,3)
                axs = fig.add_subplot(1,4,4)
                cax = fig.add_axes([.91, .1, .025, .8])
                if sigma:
                    axw.set_ylabel('sigma')
                else:
                    axw.set_ylabel('pressure')
            
                xyfactor = 1e-3 # m -> km
                     
            x = f.NCOLS + 1
            y = f.NROWS + 1
            start_south = 0
            end_south = start_south + x
            start_east = end_south
            end_east = start_east + y
            start_north = end_east
            end_north = start_north + x
            start_west = end_north
            end_west = start_west + y
            X, Y = np.meshgrid(np.arange(x + 1) * f.XCELL * xyfactor, vertcrd)
            patchess = axs.pcolor(X, Y, var[:, start_south:end_south], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
            if not map:
                if reversevert: axs.set_ylim(*axs.get_ylim()[::-1])
                axs.set_title('South')
                axs.set_xlabel('E to W km')
                axs.set_xlim(*axs.get_xlim()[::-1])
                
            X, Y = np.meshgrid(np.arange(-1, x) * f.XCELL * xyfactor, vertcrd)
            patchesn = axn.pcolor(X, Y, var[:, start_north:end_north], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
            if reversevert: axn.set_ylim(*axn.get_ylim()[::-1])
            if not map:
                axn.set_title('North')
                axn.set_xlabel('W to E km')

            if map:
                X, Y = np.meshgrid(vertcrd, np.arange(y + 1) * f.YCELL)
                patchese = axe.pcolor(X, Y, var[:, start_east:end_east].swapaxes(0,1), cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
                if reversevert: axe.set_xlim(*axe.get_xlim()[::-1])
            else:
                X, Y = np.meshgrid(np.arange(y + 1) * f.YCELL * xyfactor, vertcrd)
                patchese = axe.pcolor(X, Y, var[:, start_east:end_east], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
                if reversevert: axe.set_ylim(*axe.get_ylim()[::-1])
                axe.set_title('East')
                axe.set_xlabel('N to S km')
                axe.set_xlim(*axe.get_xlim()[::-1])
            if map:
                X, Y = np.meshgrid(vertcrd, np.arange(-1, y) * f.YCELL)
                patchesw = axw.pcolor(X, Y, var[:, start_west:end_west].swapaxes(0,1), cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
            else:
                X, Y = np.meshgrid(np.arange(-1, y) * f.YCELL * xyfactor, vertcrd)
                patchesw = axw.pcolor(X, Y, var[:, start_west:end_west], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
                if reversevert: axw.set_ylim(*axw.get_ylim()[::-1])
                axw.set_title('West')
                axw.set_xlabel('S to N km')
            if map:
                for ax in [axe, axw]:
                    ax.axis('tight', axis = 'x')
                    pl.setp( ax.xaxis.get_majorticklabels(), rotation=90 )
                for ax in [axs, axn]:
                    ax.axis('tight', axis = 'y')
            else:
                for ax in [axe, axn, axw, axs] + ([axmap] if map else []):
                    ax.axis('tight')
            
            if 'TFLAG' in f.variables.keys():
                SDATE = f.variables['TFLAG'][:][0, 0, 0]
                EDATE = f.variables['TFLAG'][:][-1, 0, 0]
                STIME = f.variables['TFLAG'][:][0, 0, 1]
                ETIME = f.variables['TFLAG'][:][-1, 0, 1]
                if SDATE == 0:
                    SDATE = 1900001
                    EDATE = 1900001
                sdate = datetime.strptime('%07d %06d' % (SDATE, STIME), '%Y%j %H%M%S')
                edate = datetime.strptime('%07d %06d' % (EDATE, ETIME), '%Y%j %H%M%S')
            elif 'tau0' in f.variables.keys():
                sdate = datetime(1985, 1, 1, 0) + timedelta(hours = f.variables['tau0'][0])
                edate = datetime(1985, 1, 1, 0) + timedelta(hours = f.variables['tau1'][-1])
            else:
                sdate = datetime(1900, 1, 1, 0)
                edate = datetime(1900, 1, 1, 0)
            try:
                title = '%s %s to %s' % (var_name, sdate.strftime('%Y-%m-%d'), edate.strftime('%Y-%m-%d'))
            except:
                title = var_name
            fig.suptitle(title.replace('O3', 'Ozone at Regional Boundaries'))
            if var.min() < vmin and var.max() > vmax:
                extend = 'both'
            elif var.min() < vmin:
                extend = 'min'
            elif var.max() > vmax:
                extend = 'max'
            else:
                extend = 'neither'
            fig.colorbar(patchesw, cax = cax, extend = extend)
            cax.set_xlabel('ppm')
            fig.savefig('%s_%s%s.%s' % (args.outpath, var_name, lstr, args.figformat))
            pl.close(fig)
    
if __name__ == '__main__':
    from PseudoNetCDF.pncparse import getparser, pncparse
    
    parser = getparser(has_ofile = True, plot_options = True, interactive = False)
    parser.add_argument("--itertime", dest = "itertime", action = "store_true", default = False,
                        help = "Create plots for each time.")

    parser.add_argument("-n", "--no-map", dest = "nomap", action = "store_true", default = False,
                        help = "Try to plot with map")

    parser.add_argument("--sigma", dest = "sigma", action = "store_true", default = False,
                        help = "Plot data on sigma coordinate.")

    parser.add_argument("--scale", dest = "scale", type = str, default = 'deciles',
                        help = "Defaults to deciles (i.e., 10 equal probability bins), but linear and log are also options.")

    parser.add_argument("--mask-zeros", dest = "maskzeros", action = "store_true", default = False,
                        help = "Defaults False.")

    parser.add_argument("--minmax", dest = "minmax", type = str, default = "None,None",
                        help = "Defaults None, None.")

    parser.add_argument("--minmaxq", dest = "minmaxq", type = str, default = '0,100',
                        help = "Defaults 0,100.")
    parser.epilog = """
Example:
    $ pncboundaries.py outputs/ts20120301.bpch.BCON.nc -v O3 test_boundary -r TSTEP,mean --states True
"""
    ifiles, args = pncparse(has_ofile = True, parser = parser)
    if args.variables is None:
        raise ValueError('User must specify variable(s) to plot:\n%s' % '\n\t'.join(ifiles[0].variables.keys()))
    
    plot(ifiles, args)
