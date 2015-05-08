#!/usr/bin/env python
# Run this script with pnc options
import os
from PseudoNetCDF.pncload import PNCConsole
from PseudoNetCDF.pncparse import pncparser
from PseudoNetCDF.coordutil import getmap, getlatbnds, getlonbnds
from warnings import warn
import numpy as np
ifiles, options = pncparser(has_ofile = True, plot_options = True)

import pylab as pl
Normalize = pl.matplotlib.colors.Normalize
LogNorm = pl.matplotlib.colors.LogNorm
BoundaryNorm = pl.matplotlib.colors.BoundaryNorm

ifile = ifiles[0]
map = getmap(ifile)
fig = pl.gcf()
ax = pl.gca()
if options.coastlines: map.drawcoastlines(ax = ax)
if options.countries: map.drawcountries(ax = ax)
if options.states: map.drawstates(ax = ax)
if options.counties: map.drawcounties(ax = ax)
nborders = len(ax.collections)

for ifile in ifiles:
    if map.projection == 'cyl':
        lat = ifile.variables['latitude']
        lon = ifile.variables['longitude']
        latb, latunit = getlatbnds(ifile)[:]
        lonb, lonunit = getlonbnds(ifile)[:]
    else:
        lat = ifile.variables['latitude']
        lon = ifile.variables['longitude']
        latb, latunit = getybnds(ifile)[:]
        lonb, lonunit = getxbnds(ifile)[:]
    
    if latb.ndim == lonb.ndim and lonb.ndim == 2:
        LON, LAT = lonb, latb
    else:
        LON, LAT = np.meshgrid(lonb, latb)
    
    variables = options.variables
    if variables is None:
        variables = [key for key, var in ifile.variables.iteritems() if len(set(['latitude', 'longitude']).intersection(getattr(var, 'coordinates', '').split())) == 2]
    if len(variables) == 0:
        raise ValueError('Unable to heuristically determin plottable variables; use -v to specify variables for plotting')
    for varkey in variables:
        ax = pl.gca()
                    
        if not options.overlay:
            del ax.collections[nborders:]
        var = ifile.variables[varkey]
        vals = var[:].squeeze()
        if options.normalize is None:
            from scipy.stats import normaltest
            vmin, vmax = vals.min(), vals.max()
            if normaltest(vals.ravel())[1] < 0.05:
                boundaries = np.logspace(np.log10(vmin), np.log10(vmax), num = 10)
                warn('Selected LogScale for %s' % varkey)
            else:
                boundaries = np.linspace(vmin, vmax, num = 10)
            norm = BoundaryNorm(boundaries, ncolors = 256)
        else:
            norm = eval(options.normalize)
        varunit = getattr(var, 'units', 'unknown').strip()
        print varkey,
        if vals.ndim == 1:
            patches = map.scatter(lon[:], lat[:], c = vals, s = 24, norm = norm, ax = ax)
        else:
            patches = map.pcolor(LON, LAT, vals, norm = norm, ax = ax)
        if lonunit == 'x (LCC m)':
            ax.xaxis.get_major_formatter().set_scientific(True)
            ax.xaxis.get_major_formatter().set_powerlimits((-3, 3))
        if latunit == 'y (LCC m)':
            ax.yaxis.get_major_formatter().set_scientific(True)
            ax.yaxis.get_major_formatter().set_powerlimits((-3, 3))
        ax.set_xlabel(lonunit)
        ax.set_ylabel(latunit)
        height = LAT.max() - LAT.min()
        width = LON.max() - LON.min()
        if width > height:
            orientation = 'horizontal'
        else:
            orientation = 'vertical'
        try:
            cax = cbar.ax
        except:
            cax = None
        cbar = pl.gcf().colorbar(patches, orientation = orientation, cax = cax)
        del cbar.ax.texts[:]
        cbar.set_label(varkey + ' (' + varunit + ')')
        if orientation == 'vertical':
            cbar.ax.text(.5, 1, '%.2g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
            cbar.ax.text(.5, 0, '%.2g ' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
        else:
            cbar.ax.text(1, .5, ' %.2g' % var[:].max(), verticalalignment = 'center', horizontalalignment = 'left')
            cbar.ax.text(0, .5, '%.2g ' % var[:].min(), verticalalignment = 'center', horizontalalignment = 'right')
        try:
            cbar.formatter.set_scientific(True)
            cbar.formatter.set_powerlimits((-3, 3))
        except:
            pass
        cbar.update_ticks()
        fmt = 'png'
        figpath = os.path.join(options.outpath + '_map_' + varkey + '.' + fmt)
        if options.interactive:
            csl = PNCConsole(locals = globals())
            csl.interact()
            
        pl.savefig(figpath)
        print 'Saved fig', figpath
        
