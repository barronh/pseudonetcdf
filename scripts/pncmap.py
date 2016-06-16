#!/usr/bin/env python
from __future__ import print_function
# Run this script with pnc options
import sys
import os

from PseudoNetCDF.pncload import PNCConsole
from PseudoNetCDF.pncparse import pncparse, getparser
from PseudoNetCDF.coordutil import getmap, getlatbnds, getlonbnds, getybnds, getxbnds
import warnings
warn=warnings.warn
import numpy as np
from PseudoNetCDF.plotutil import *
def prepmap(ifiles, args):
    ifile = ifiles[0]
    ax = pl.gca()
    map = getmap(ifile, resolution = args.resolution)
    if args.coastlines: map.drawcoastlines(ax = ax)
    if args.countries: map.drawcountries(ax = ax)
    if args.states: map.drawstates(ax = ax)
    if args.counties: map.drawcounties(ax = ax)
    for si, shapefile in enumerate(args.shapefiles):
        shapeopts = shapefile.split(',')
        shapepath = shapeopts[0]
        shapeoptdict = eval('dict(' + ','.join(shapeopts[1:]) + ')')
        shapename = os.path.basename(shapepath)[:-3] + str(si)
        map.readshapefile(shapepath, shapename, ax = ax, **shapeoptdict)
    args.map = map

def makemap(ifiles, args):
    fig = pl.gcf()
    if len(args.figure_keywords) > 0:
        plt.setp(fig, **args.figure_keywords)
    
    ax = pl.gca()
    if len(args.axes_keywords) > 0:
        plt.setp(ax, **args.axes_keywords)
    
    map = args.map
    nborders = len(ax.collections)
    for fi, ifile in enumerate(ifiles):
        if map.projection in ('lcc', 'merc'):
            lat = ifile.variables['latitude']
            lon = ifile.variables['longitude']
            latb, latunit = getybnds(ifile)[:]
            lonb, lonunit = getxbnds(ifile)[:]
        else:
            lat = ifile.variables['latitude']
            lon = ifile.variables['longitude']
            latb, latunit = getlatbnds(ifile)[:]
            lonb, lonunit = getlonbnds(ifile)[:]
        
        if latb.ndim == lonb.ndim and lonb.ndim == 2:
            LON, LAT = lonb, latb
        else:
            LON, LAT = np.meshgrid(lonb.view(np.ndarray), latb.view(np.ndarray))
    
        variables = args.variables
        if variables is None:
            variables = [key for key, var in ifile.variables.items() if len(set(['latitude', 'longitude']).intersection(getattr(var, 'coordinates', '').split())) == 2]
        if len(variables) == 0:
            raise ValueError('Unable to heuristically determin plottable variables; use -v to specify variables for plotting')
        for varkey in variables:
            ax = pl.gca()
                    
            if not args.overlay:
                del ax.collections[nborders:]
            var = ifile.variables[varkey]
            if args.squeeze:
                vals = var[:].squeeze()
            else:
                vals = var[:]
            vmin, vmax = vals.min(), vals.max()
            if args.normalize is None:
                from scipy.stats import normaltest
                if normaltest(vals.ravel())[1] < 0.001:
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

            if not norm.vmin is None:
                vmin = norm.vmin
            if not norm.vmax is None:
                vmax = norm.vmax
            varunit = getattr(var, 'units', 'unknown').strip()
            if args.verbose > 0: print(varkey, sep = '')
            if vals.ndim == 1:
                patches = map.scatter(lon[:], lat[:], c = vals, edgecolors = 'none', s = 24, norm = norm, ax = ax, zorder = 2)
            else:
                patches = map.pcolor(LON, LAT, vals, norm = norm, ax = ax)
            if lonunit == 'x (m)':
                ax.xaxis.get_major_formatter().set_scientific(True)
                ax.xaxis.get_major_formatter().set_powerlimits((-3, 3))
            if latunit == 'y (m)':
                ax.yaxis.get_major_formatter().set_scientific(True)
                ax.yaxis.get_major_formatter().set_powerlimits((-3, 3))
            ax.set_xlabel(lonunit)
            ax.set_ylabel(latunit)
            height = np.abs(np.diff(ax.get_ylim()))
            width = np.abs(np.diff(ax.get_xlim()))
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
            cbar = pl.gcf().colorbar(patches, orientation = orientation, cax = cax, extend = extend, format = formatter, spacing = 'proportional')
            del cbar.ax.texts[:]
            cbar.set_label(varkey + ' (' + varunit + '; min=%.3g; max=%.3g)' % (var[:].min(), var[:].max()))
 #           if orientation == 'vertical':
 #               cbar.ax.text(.5, 1.05, '%.3g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
#                cbar.ax.text(.5, -.06, '%.3g ' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
#            else:
#                cbar.ax.text(1.05, .5, ' %.3g' % var[:].max(), verticalalignment = 'center', horizontalalignment = 'left')
#                cbar.ax.text(-.06, .5, '%.3g ' % var[:].min(), verticalalignment = 'center', horizontalalignment = 'right')
            cbar.update_ticks()
            fmt = args.figformat
            outpath = args.outpath
            if len(ifiles) > 1:
                lstr = str(fi).rjust(len(str(len(ifiles))), '0')
                if args.verbose > 0: print('adding numeric suffix for file', lstr)
            else:
                lstr = ''
                
            figpath = os.path.join(outpath + varkey + lstr + '.' + fmt)
            if args.interactive:
                csl = PNCConsole(locals = globals())
                csl.interact()
            for cmd in args.plotcommands:
                exec(cmd)
            pl.savefig(figpath)
            if args.verbose > 0: print('Saved fig', figpath)
        
if __name__ == '__main__':
    parser = getparser(has_ofile = True, map_options = True, plot_options = True, interactive = True)
    parser.add_argument('--no-squeeze', dest = 'squeeze', default = True, action = 'store_false', help = 'Squeeze automatically removes singleton dimensions; disabling requires user to remove singleton dimensions with --remove-singleton option')
    parser.add_argument("--iter", dest = "iter", action = 'append', default = [], help = "Create plots for each element of specified dimension (e.g., --iter=time).")
    ifiles, args = pncparse(has_ofile = True, plot_options = True, interactive = True, parser = parser)
    prepmap(ifiles, args)
    if args.iter != []:
        ifile, = ifiles
        ifiles = []
        for dimk in args.iter:
            ifiles += [slice_dim(getvarpnc(ifile, None), '%s,%d' % (dimk,i)) for i in range(len(ifile.dimensions[dimk]))]
    makemap(ifiles, args)