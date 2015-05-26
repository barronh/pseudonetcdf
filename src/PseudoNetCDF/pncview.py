__all__ = ['mapplot', 'profile', 'tileplot', 'presslon', 'presslat', 'timeseries', 'OptionDict']

from warnings import warn
from Tkinter import Checkbutton, Frame, Label, Scrollbar, Listbox, Button, IntVar, Tk, VERTICAL, EXTENDED, END, N, S, SINGLE, Entry, StringVar, Text, DISABLED, PhotoImage, LEFT, E, W
import os
from types import MethodType
import pylab as pl

import numpy as np
_pre_code = 'pl.close(); pl.rcParams["image.cmap"] = "jet"'
_before_code = 'pl.figure(); pl.interactive(True); '
_after_code = 'ax.set_title(varkey); pl.draw()'
_post_code = ''
_coastlines_opt = True
_countries_opt = True
_states_opt = True
_counties_opt = False

from PseudoNetCDF.coordutil import *
class OptionDict(dict):
    def __init__(self, *args, **kwds):
        """
        coastlines = %(_coastlines_opt)s
        countries = %(_countries_opt)s
        states = %(_states_opt)s
        counties = %(_counties_opt)s
        pre_txt = Code to run before any plots (%(_pre_code)s)
        before_txt = Code to run before each variable-plot (%(_before_code)s)
        after_txt = Code to run after each plot (%(_after_code)s)
        post_txt = Code to run after all plots (%(_post_code)s)
        """ % globals()
        dict.__init__(self, coastlines = _coastlines_opt, countries = _countries_opt, states = _states_opt, counties = _counties_opt, pre_txt = _pre_code, before_txt = _before_code, after_txt = _after_code, post_txt = _post_code)
        dict.__init__(self, *args, **kwds)
        if 'outpath' not in self:
            raise KeyError('outpath is a required option')
    def __getattr__(self, key):
        return self.get(key, False)

defaultoption = OptionDict(outpath = 'pncview')

class TkApp:
    def __init__(self, ncffile, options):
        master = self.root = Tk()
        self.ncffile = ncffile
        self.options = options
        self.plotted_variables = set()
        frame = Frame(master)
        frame.grid(row = 0)
        codeframe = Frame(master)
        codeframe.grid(row = 1)
        metaframe = Frame(master)
        metaframe.grid(row = 2)
        goframe = Frame(frame)
        goframe.grid(column = 3, row = 1)

        var_label = Label(frame, text = 'Select Variable')
        var_label.grid(column = 0, row = 0)
        var_scrollbar = Scrollbar(frame, orient = VERTICAL)
        var_scrollbar.grid(column = 1, row = 1, sticky = N + S)
        self.var = Listbox(frame, selectmode = EXTENDED, exportselection = 0, yscrollcommand = var_scrollbar.set)
        self.var.grid(column = 0, row = 1)
        var_scrollbar.config(command = self.var.yview)

        what_to_do = Label(frame, text = 'Execute')
        what_to_do.grid(column = 2, row = 0)
        self.method_list = Listbox(frame, selectmode = SINGLE, exportselection = 0)
        self.method_list.grid(column = 2, row = 1)
        self.pre_txt = StringVar()
        pre_label = Label(codeframe, text = 'Before any figures, execute code')
        self.pre_txt.set(_pre_code)
        pre_label.grid(row = 2, sticky = 'W')
        self.pre = Entry(codeframe, width = 120, textvariable = self.pre_txt)
        self.pre.grid(row =3, sticky = 'E')

        self.before_txt = StringVar()
        self.before_txt.set(_before_code)
        before_label = Label(codeframe, text = 'Before each figure, execute code')
        before_label.grid(row = 4, sticky = 'W')
        self.before = Entry(codeframe, width = 120, textvariable = self.before_txt)
        self.before.grid(row =5, sticky = 'E')

        self.after_txt = StringVar()
        self.after_txt.set(_after_code)
        after_label = Label(codeframe, text = 'After each figure, execute code')
        after_label.grid(row = 6, sticky = 'W')
        self.after = Entry(codeframe, width = 120, textvariable = self.after_txt)
        self.after.grid(row =7, sticky = 'E')

        self.post_txt = StringVar()
        self.post_txt.set(_post_code) 
        post_label = Label(codeframe, text = 'After all figures, execute code')
        post_label.grid(row = 8, sticky = 'W')
        self.post = Entry(codeframe, width = 120, textvariable = self.post_txt)
        self.post.grid(row = 9, sticky = 'E')

        options_label = Label(goframe, text = 'Options:')
        options_label.grid(column = 0, row = 1, sticky = 'W')
        self.logscale = IntVar()
        self.logscale.set(0)
        c = Checkbutton(goframe, text = "log-scale?", variable = self.logscale)
        c.grid(column = 0, row = 2, sticky = 'W')

        self.coastlines = IntVar()
        self.coastlines.set(_coastlines_opt)
        coastlines = Checkbutton(goframe, text = "coastlines?", variable = self.coastlines, justify = LEFT)
        coastlines.grid(column = 0, row = 3, sticky = 'W')

        self.countries = IntVar()
        self.countries.set(_countries_opt)
        countries = Checkbutton(goframe, text = "countries?", variable = self.countries, justify = LEFT)
        countries.grid(column = 0, row = 4, sticky = 'W')

        self.states = IntVar()
        self.states.set(_states_opt)
        states = Checkbutton(goframe, text = "states?", variable = self.states, justify = LEFT)
        states.grid(column = 0, row = 5, sticky = 'W')

        self.counties = IntVar()
        self.counties.set(_counties_opt)
        counties = Checkbutton(goframe, text = "counties?", variable = self.counties, justify = LEFT)
        counties.grid(column = 0, row = 6, sticky = 'W')

        self.execute_button = Button(goframe, text = "Make Figure", command = self.execute)
        self.execute_button.grid(row = 0, column = 0, sticky = 'W')
                
        self.methods = ['mapplot', 'presslat', 'presslon', 'time-lat', 'profile', 'timeseries', 'pressx', 'tileplot', 'plot']
        method_labels= ['lat-lon', 'press-lat', 'press-lon', 'time-lat', 'Vertical Profile', 'Time Series', 'press-? (2-D)', 'Tile Plot (2-D)', 'Plot (1-D)']
        for method in method_labels:
            self.method_list.insert(END, method)
            
        var_keys = [k for k, v in self.ncffile.variables.iteritems() if k not in ('time', 'latitude', 'longitude', 'latitude_bounds', 'longitude_bounds', 'time_bounds', 'tau0', 'tau1', 'TFLAG')]
        var_keys.sort()
        self.vars = []
        for spc in var_keys:
            self.var.insert(END, spc)
            self.vars.append(spc)

        meta_label = Label(metaframe, text = 'Common Data Language Header:')
        meta_label.grid(column = 0, row = 0, sticky = 'W')
        meta_scrollbar = Scrollbar(metaframe, orient = VERTICAL)
        meta_scrollbar.grid(column = 1, row = 1, sticky = N + S)
        self.meta = Text(metaframe, height=10, width=118, bg='white', relief='flat', yscrollcommand = meta_scrollbar.set)
        self.meta.grid(column = 0, row = 1, sticky = 'W')
        from PseudoNetCDF.pncdump import pncdump
        from StringIO import StringIO
        pdump = StringIO("")
        try:
            name = ', '.join(options.ifile)
        except:
            name = 'ifile'
        pncdump(self.ncffile, header = True, outfile = pdump, )
        pdump.seek(0, 0)
        self.meta.insert(END, pdump.read())
        self.meta.config(state=DISABLED)
        help = Button(goframe, text = 'Help', command = self.help)
        help.grid(column = 0, row = 7)
        quit = Button(goframe, text = 'Quit', command = self.quit)
        quit.grid(column = 0, row = 8)
        master.mainloop()
   
    def help(self):
        print "pl is pylab: details at matplotlib;"
        
    def quit(self):
        self.root.destroy()

    def _get_var(self, list):
        items = list.curselection()
        try: items = map(int, items)
        except: pass
        items = [self.vars[i] for i in items]
        return items
        
    def get_var(self):
        return self._get_var(self.var)
        
    def get_methods(self):
        items = self.method_list.curselection()
        try: items = map(int, items)
        except: pass
        items = [self.methods[i] for i in items]
        return items
    def execute(self):
        os.system('clear')
        vars = self.get_var()
        self.plotted_variables = self.plotted_variables.union(vars)
        methods, = self.get_methods()
        self.options.logscale = bool(self.logscale.get())
        self.options.coastlines = bool(self.coastlines.get())
        self.options.countries = bool(self.countries.get())
        self.options.states = bool(self.states.get())
        self.options.counties = bool(self.counties.get())
        self.options.pre_txt = self.pre_txt.get()
        self.options.before_txt = self.before_txt.get()
        self.options.after_txt = self.after_txt.get()
        self.options.post_txt = self.post_txt.get()
        plotwithopts(self.ncffile, methods, vars, self.options)

def plotwithopts(ifile, method, vars, options = defaultoption):
    from PseudoNetCDF.sci_var import getvarpnc
    from PseudoNetCDF.pncgen import pncgen
    exec(options.pre_txt)
    for varkey in vars:
        figpath = eval(method)(ifile = ifile, varkey = varkey, options = options, before = options.before_txt, after = options.after_txt)
        pncgen(getvarpnc(ifile, list(vars) + ['TFLAG', 'time', 'latitude', 'longitude', 'latitude_bounds', 'longitude_bounds']), figpath + '.nc', verbose = False)
    exec(options.post_txt)

def StartTk(ncffile, options):
    TkApp(ncffile, options)

def pncview(ifile, options):
    # add a gui for plotting
    return StartTk(ifile, options)

def gettime(ifile):
    from PseudoNetCDF import PseudoNetCDFVariable
    from pylab import date2num
    from datetime import datetime
    if 'time' in ifile.variables:
        time = ifile.variables['time']
        unit = time.units
        tunit, datestr = unit.split(' since ')
        sdate = datetime.strptime(datestr.replace(' UTC', ''), '%Y-%m-%d %H:%M:%S')
        time = np.array([sdate + timedelta(**{tunit: t}) for t in time[:]])
        unit = 'time'
    elif 'TFLAG' in ifile.variables:
        x = ifile.variables['TFLAG'][:].copy()
        for axi in range(1, x.ndim - 1):
            x = x[:, 0]
        time = np.array([(datetime.strptime('%07d %06d' % (d, t), '%Y%j %H%M%S')) for d, t in x])
        unit = 'time'
    elif 'time' in ifile.dimensions:
        time = np.arange(len(ifile.dimensions['time']))
        unit = 'steps'
    elif 'TSTEP' in ifile.dimensions:
        time = np.arange(1, len(ifile.dimensions['TSTEP']) + 1)
        unit = 'steps'
    else:
        raise KeyError('No time found')
    return PseudoNetCDFVariable(None, 'time', 'f', ('time',), values = time[:], units = unit)
    
def timeseries(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    outpath = getattr(options, 'outpath', '.')
    time = gettime(ifile)
    var = ifile.variables[varkey]
    dims = [(k, l) for l, k in zip(var[:].shape, var.dimensions) if l > 1]
    if len(dims) > 1:
        raise ValueError('Time series can have 1 non-unity dimensions; got %d - %s' % (len(dims), str(dims)))
    exec(before)
    ax = pl.gca()
    print varkey,
    if options.logscale:
        ax.set_yscale('log')
        
    ax.plot_date(time[:].squeeze(), var[:].squeeze())
    ax.set_xlabel(time.units.strip())
    ax.set_ylabel(getattr(var, 'standard_name', varkey).strip() + ' ' + var.units.strip())
    fmt = 'png'
    figpath = os.path.join(outpath + '_ts_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def plot(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    outpath = getattr(options, 'outpath', '.')
    var = ifile.variables[varkey]
    dims = [(k, l) for l, k in zip(var[:].shape, var.dimensions) if l > 1]
    if len(dims) > 1:
        raise ValueError('Plots can have only 1 non-unity dimensions; got %d - %s' % (len(dims), str(dims)))
    exec(before)
    ax = pl.gca()
    print varkey,
    if options.logscale:
        ax.set_yscale('log')
        
    ax.plot(var[:].squeeze())
    ax.set_xlabel('unknown')
    ax.set_ylabel(getattr(var, 'standard_name', varkey).strip() + ' ' + var.units.strip())
    fmt = 'png'
    figpath = os.path.join(outpath + '_1d_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def pressx(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    from matplotlib.colors import Normalize, LogNorm
    outpath = getattr(options, 'outpath', '.')
    vert = getpresbnds(ifile)
    var = ifile.variables[varkey]
    dims = [(k, l) for l, k in zip(var[:].shape, var.dimensions) if l > 1]
    if len(dims) > 2:
        raise ValueError('Press-x can have 2 non-unity dimensions; got %d - %s' % (len(dims), str(dims)))
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    exec(before)
    ax = pl.gca()
    print varkey,
    vals = var[:].squeeze()
    x = np.arange(vals.shape[1])
    patches = ax.pcolor(x, vert, vals, norm = norm)
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    pl.colorbar(patches)
    ax.set_ylim(vert.max(), vert.min())
    ax.set_xlim(x.min(), x.max())
    fmt = 'png'
    figpath = os.path.join(outpath + '_PRESX_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def presslat(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    from matplotlib.colors import Normalize, LogNorm
    outpath = getattr(options, 'outpath', '.')
    vert = getpresbnds(ifile)
    lat, latunit = getlatbnds(ifile)
    lat = np.append(lat.squeeze()[..., :2].mean(1), lat.squeeze()[-1, 2:].mean(0))
    var = ifile.variables[varkey]
    dims = [(k, l) for l, k in zip(var[:].shape, var.dimensions) if l > 1]
    if len(dims) > 2:
        raise ValueError('Press-lat can have 2 non-unity dimensions; got %d - %s' % (len(dims), str(dims)))
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    exec(before)
    ax = pl.gca()
    print varkey,
    patches = ax.pcolor(lat, vert, var[:].squeeze(), norm = norm)
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    cbar = pl.colorbar(patches)
    vunit = getattr(var, 'units', 'unknown').strip()
    cbar.set_label(varkey + ' (' + vunit + ')')
    cbar.ax.text(.5, 1, '%.2g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
    cbar.ax.text(.5, 0, '%.2g' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
    ax.set_ylim(vert.max(), vert.min())
    ax.set_xlim(lat.min(), lat.max())
    fmt = 'png'
    figpath = os.path.join(outpath + '_PRESSLAT_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def presslon(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    from matplotlib.colors import Normalize, LogNorm
    outpath = getattr(options, 'outpath', '.')
    vert = getpresbnds(ifile)
    lon, lonunit = getlonbnds(ifile)
    lon = np.append(lon.squeeze()[..., [0, 3]].mean(1), lon.squeeze()[-1, [1, 2]].mean(0))
    var = ifile.variables[varkey]
    dims = [(k, l) for l, k in zip(var[:].shape, var.dimensions) if l > 1]
    if len(dims) > 2:
        raise ValueError('Press-lon plots can have 2 non-unity dimensions; got %d - %s' % (len(dims), str(dims)))
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    exec(before)
    ax = pl.gca()
    print varkey,
    patches = ax.pcolor(lon, vert, var[:].squeeze(), norm = norm)
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    cbar = pl.colorbar(patches)
    vunit = getattr(var, 'units', 'unknown').strip()
    cbar.set_label(varkey + ' (' + vunit + ')')
    cbar.ax.text(.5, 1, '%.2g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
    cbar.ax.text(.5, 0, '%.2g' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')

    ax.set_ylim(vert.max(), vert.min())
    ax.set_xlim(lon.min(), lon.max())
    fmt = 'png'
    figpath = os.path.join(outpath + '_PRESLON_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def tileplot(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    from matplotlib.colors import Normalize, LogNorm
    outpath = getattr(options, 'outpath', '.')
    var = ifile.variables[varkey]
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    exec(before)
    ax = pl.gca()
    print varkey,
    dims = [(k, l) for l, k in zip(var[:].shape, var.dimensions) if l > 1]
    if len(dims) > 2:
        raise ValueError('Tile plots can have 2 non-unity dimensions; got %d - %s' % (len(dims), str(dims)))
    patches = ax.pcolor(var[:].squeeze(), norm = norm)
    ax.set_xlim(0, var[:].squeeze().shape[1])
    ax.set_ylim(0, var[:].squeeze().shape[0])
    ax.set_xlabel(dims[1][0])
    ax.set_ylabel(dims[0][0])
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    cbar = pl.colorbar(patches)
    vunit = getattr(var, 'units', 'unknown').strip()
    cbar.set_label(varkey + ' (' + vunit + ')')
    cbar.ax.text(.5, 1, '%.2g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
    cbar.ax.text(.5, 0, '%.2g' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
    fmt = 'png'
    figpath = os.path.join(outpath + '_2D_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def minmaxmean(ax, vals, vertcrd, **kwds):
    minval = vals.min(1)
    meanval = vals.mean(1)
    maxval = vals.max(1)
    kwds.setdefault('edgecolor', 'k')
    kwds.setdefault('facecolor', 'grey')
    kwds.setdefault('alpha', 1)
    kwds.setdefault('ls', 'solid')
    linekwds = kwds.copy()
    linekwds['color'] = linekwds.pop('edgecolor')
    linekwds.pop('facecolor')
    linekwds.pop('alpha')
    linekwds.pop('ls')
    fillkwds = kwds.copy()
    
    line, = ax.plot(meanval, vertcrd, **linekwds)

    x = np.ma.concatenate([minval[:vertcrd.size], maxval[:vertcrd.size][::-1]])
    y = np.ma.concatenate([vertcrd[:], vertcrd[::-1]])
    mask = x.mask | y.mask
    x = np.ma.masked_where(mask, x).compressed()
    y = np.ma.masked_where(mask, y).compressed()
    range, = ax.fill(x, y, **fillkwds)
    return line, range

def profile(ifile, varkey, options, before = '', after = ''):
    import pylab as pl
    print varkey,
    outpath = getattr(options, 'outpath', '.')
    try:
        vert = getpresmid(ifile)
        vunit = 'Pa'
    except:
        vert = getsigmamid(ifile)
        vunit = r'\sigma'
    var = ifile.variables[varkey]
    
    dims = list(var.dimensions)
    for knownvert in ['layer', 'LAY'] + ['layer%d' % i for i in range(72)]:
        if knownvert in dims:
            vidx = dims.index(knownvert)
            break
    else:
        raise KeyError("No known vertical coordinate; got %s" % str(dims))
    vert = vert[:var[:].shape[vidx]]
    units = var.units.strip()
    vals = np.rollaxis(var[:], vidx, start = 0).view(np.ma.MaskedArray).reshape(vert.size, -1)
    ax = pl.gca()
    minmaxmean(ax, vals, vert)
    ax.set_xlabel(varkey + ' ('+units+')')
    ax.set_ylabel(vunit)
    ax.set_ylim(vert.max(), vert.min())
    if options.logscale: ax.set_xscale('log')
    fmt = 'png'
    figpath = os.path.join(outpath + '_profile_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath
    
def mapplot(ifile, varkey, options, before = '', after = ''):
    """
    ifile - a pseudonetcdf file
    varkey - the variable to plot
    options - 
    """
    import pylab as pl
    from matplotlib.colors import Normalize, LogNorm
    outpath = getattr(options, 'outpath', '.')
    map = getmap(ifile)
    if map.projection == 'cyl':
        latb, latunit = getlatbnds(ifile)[:]
        lonb, lonunit = getlonbnds(ifile)[:]
    else:
        latb, latunit = getybnds(ifile)[:]
        lonb, lonunit = getxbnds(ifile)[:]
    if latb.ndim == lonb.ndim and lonb.ndim == 2:
        LON, LAT = lonb, latb
    else:
        LON, LAT = np.meshgrid(lonb, latb)

    ax = pl.gca()
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    var = ifile.variables[varkey]
    exec(before)
    ax = pl.gca()
    vunit = getattr(var, 'units', 'unknown').strip()
    print varkey,
    try:
        if options.coastlines: map.drawcoastlines(ax = ax)
        if options.countries: map.drawcountries(ax = ax)
        if options.states: map.drawstates(ax = ax)
        if options.counties: map.drawcounties(ax = ax)
    except:
        print 'nomap'
        pass
    patches = map.pcolor(LON, LAT, var[:].squeeze(), norm = norm, ax = ax)
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
    cbar = pl.gcf().colorbar(patches, orientation = orientation)
    cbar.set_label(varkey + ' (' + vunit + ')')
    if orientation == 'vertical':
        cbar.ax.text(.5, 1, '%.2g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
        cbar.ax.text(.5, 0, '%.2g' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
    else:
        cbar.ax.text(1, .5, '%.2g' % var[:].max(), verticalalignment = 'center', horizontalalignment = 'left')
        cbar.ax.text(0, .5, '%.2g' % var[:].min(), verticalalignment = 'center', horizontalalignment = 'right')
    try:
        cbar.formatter.set_scientific(True)
        cbar.formatter.set_powerlimits((-3, 3))
    except:
        pass
    cbar.update_ticks()
    fmt = 'png'
    figpath = os.path.join(outpath + '_map_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath
            
def main():
    import pylab as pl
    from pncparse import pncparse
    ifiles, options = pncparse(has_ofile = True, plot_options = True, interactive = False)
    if len(ifiles) != 1:
        raise IOError('pncview can operate on only 1 file; user requested %d' % len(ifiles))
    ifile, = ifiles
    pl.interactive(True)
    for method_vars in options.plotcommands:
        pieces = method_vars.split(',')
        plotargs = [p for p in pieces if '=' not in p]
        plotkwds = [p for p in pieces if '=' in p]
        method, = plotargs[:1]
        vars = plotargs[1:]
        plotoptions = eval('OptionDict(outpath="%s",%s)' % (options.outpath, ','.join(plotkwds)))
        print plotoptions.logscale
        plotwithopts(ifile, method, vars, plotoptions)
    pl.interactive(False)
    if len(options.plotcommands) == 0:
        pncview(ifile, options)

if __name__ == '__main__':
    main()
