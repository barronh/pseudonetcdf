#!/usr/bin/env python
from __future__ import print_function
# Run this script with pnc options
import sys
import os

from PseudoNetCDF.pncload import PNCConsole
from PseudoNetCDF.coordutil import getmap, getlatbnds, getlonbnds, getybnds, getxbnds
from PseudoNetCDF.sci_var import getvarpnc, slice_dim

import warnings
warn=warnings.warn
import numpy as np
from PseudoNetCDF.plotutil import *

def makemaps(args):
    ifiles = args.ifiles
    ifile = ifiles[0]
    if args.iter != []:
        ifile, = ifiles
        ifiles = []
        for dimk in args.iter:
            ifiles += [slice_dim(getvarpnc(ifile, None), '%s,%d' % (dimk,i)) for i in range(len(ifile.dimensions[dimk]))]
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
                notmasked = ~(np.ma.getmaskarray(lon[:]) | np.ma.getmaskarray(lat[:]) | np.ma.getmaskarray(vals[:]))
                scatlon = lon[:][notmasked]
                scatlat = lat[:][notmasked]
                scatvals = vals[:][notmasked]
                patches = map.scatter(scatlon[:], scatlat[:], c = scatvals, edgecolors = 'none', s = 24, norm = norm, ax = ax, zorder = 2)
            else:
                if vals.ndim != 2:
                    warn('Maps require 2-d data; values right now %s %s' % (str(vals.shape), str(dict(zip(var.dimensions, var.shape)))))
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
