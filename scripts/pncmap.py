#!/usr/bin/env python
# Run this script with pnc options

import sys
import os

from PseudoNetCDF.pncload import PNCConsole
from PseudoNetCDF.pncparse import pncparse
from PseudoNetCDF.coordutil import getmap, getlatbnds, getlonbnds, getybnds, getxbnds
import warnings
warn=warnings.warn
import numpy as np
from PseudoNetCDF.plotutil import *
def prepmap(ifiles, options):
    ifile = ifiles[0]
    ax = pl.gca()
    map = getmap(ifile, resolution = 'c')
    if options.coastlines: map.drawcoastlines(ax = ax)
    if options.countries: map.drawcountries(ax = ax)
    if options.states: map.drawstates(ax = ax)
    if options.counties: map.drawcounties(ax = ax)
    for si, shapefile in enumerate(options.shapefiles):
        shapename = os.path.basename(shapefile)[:-3] + str(si)
        map.readshapefile(shapefile, shapename, ax = ax, linewidth = pl.rcParams['lines.linewidth'])
    options.map = map

def makemap(ifile, options):
    fig = pl.gcf()
    ax = pl.gca()
    map = options.map
    nborders = len(ax.collections)
    for fi, ifile in enumerate(ifiles):
        if map.projection == 'lcc':
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
            vmin, vmax = vals.min(), vals.max()
            if options.normalize is None:
                from scipy.stats import normaltest
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
                norm = eval(options.normalize)
                formatter = None
            if not options.colorbarformatter is None:
                try:
                    formatter = eval(options.colorbarformatter)
                except:
                    formatter = options.colorbarformatter

                
            vmin, vmax = norm.vmin, norm.vmax
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
            cbar = pl.gcf().colorbar(patches, orientation = orientation, cax = cax, extend = extend, format = formatter)
            del cbar.ax.texts[:]
            cbar.set_label(varkey + ' (' + varunit + '; min=%.3g; max=%.3g)' % (var[:].min(), var[:].max()))
 #           if orientation == 'vertical':
 #               cbar.ax.text(.5, 1.05, '%.3g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
#                cbar.ax.text(.5, -.06, '%.3g ' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
#            else:
#                cbar.ax.text(1.05, .5, ' %.3g' % var[:].max(), verticalalignment = 'center', horizontalalignment = 'left')
#                cbar.ax.text(-.06, .5, '%.3g ' % var[:].min(), verticalalignment = 'center', horizontalalignment = 'right')
            cbar.update_ticks()
            fmt = 'png'
            outpath = options.outpath
            if len(ifiles) > 1:
                
                outpath += ('_%%0%dd' % len(str(len(ifiles)))) % fi
            figpath = os.path.join(outpath + varkey + '.' + fmt)
            if options.interactive:
                csl = PNCConsole(locals = globals())
                csl.interact()
            
            pl.savefig(figpath)
            print 'Saved fig', figpath
        
if __name__ == '__main__':
    ifiles, options = pncparse(has_ofile = True, plot_options = True, interactive = True)
    prepmap(ifiles, options)
    makemap(ifiles, options)
