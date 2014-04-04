__all__ = ['mapplot', 'profile', 'tileplot', 'presslon', 'presslat', 'timeseries', 'OptionDict']

from matplotlib import use; use('TkAgg')
from Tkinter import Checkbutton, Frame, Label, Scrollbar, Listbox, Button, IntVar, Tk, VERTICAL, EXTENDED, END, N, S, SINGLE, Entry, StringVar, Text, DISABLED, PhotoImage, LEFT, E, W
import os
from types import MethodType
import pylab as pl
from matplotlib.colors import Normalize, LogNorm
import numpy as np

class OptionDict(dict):
    def __init__(self, *args, **kwds):
        dict.__init__(self, *args, **kwds)
        if 'outpath' not in self:
            raise KeyError('outpath is a required option')
    def __getattr__(self, key):
        return self.get(key, False)

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
        self.pre_txt.set('pl.figure(); pl.rcParams["image.cmap"] = "jet"')
        pre_label.grid(row = 2, sticky = 'W')
        self.pre = Entry(codeframe, width = 120, textvariable = self.pre_txt)
        self.pre.grid(row =3, sticky = 'E')

        self.before_txt = StringVar()
        self.before_txt.set('pl.clf();')
        before_label = Label(codeframe, text = 'Before each figure, execute code')
        before_label.grid(row = 4, sticky = 'W')
        self.before = Entry(codeframe, width = 120, textvariable = self.before_txt)
        self.before.grid(row =5, sticky = 'E')

        self.after_txt = StringVar()
        self.after_txt.set('ax.set_title(varkey); pl.show()')
        after_label = Label(codeframe, text = 'After each figure, execute code')
        after_label.grid(row = 6, sticky = 'W')
        self.after = Entry(codeframe, width = 120, textvariable = self.after_txt)
        self.after.grid(row =7, sticky = 'E')

        self.post_txt = StringVar()
        self.post_txt.set('pl.close()') 
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
        self.coastlines.set(0)
        coastlines = Checkbutton(goframe, text = "coastlines?", variable = self.coastlines, justify = LEFT)
        coastlines.grid(column = 0, row = 3, sticky = 'W')

        self.countries = IntVar()
        self.countries.set(0)
        countries = Checkbutton(goframe, text = "countries?", variable = self.countries, justify = LEFT)
        countries.grid(column = 0, row = 4, sticky = 'W')

        self.states = IntVar()
        self.states.set(0)
        states = Checkbutton(goframe, text = "states?", variable = self.states, justify = LEFT)
        states.grid(column = 0, row = 5, sticky = 'W')

        self.counties = IntVar()
        self.counties.set(0)
        counties = Checkbutton(goframe, text = "counties?", variable = self.counties, justify = LEFT)
        counties.grid(column = 0, row = 6, sticky = 'W')

        self.execute_button = Button(goframe, text = "Make Figure", command = self.execute)
        self.execute_button.grid(row = 0, column = 0, sticky = 'W')
                
        self.methods = ['mapplot', 'presslat', 'presslon', 'time-lat', 'profile', 'timeseries', 'tileplot', 'plot']
        method_labels= ['lat-lon', 'press-lat', 'press-lon', 'time-lat', 'Vertical Profile', 'Time Series', 'Tile Plot (2-D)', 'Plot (1-D)']
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
        pncdump(self.ncffile, header = True, outfile = pdump, name = ', '.join(options.ifile))
        pdump.seek(0, 0)
        self.meta.insert(END, pdump.read())
        self.meta.config(state=DISABLED)
        quit = Button(goframe, text = 'Quit', command = self.quit)
        quit.grid(column = 0, row = 7)
        master.mainloop()
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
        from PseudoNetCDF.pncgen import pncgen
        from PseudoNetCDF.sci_var import getvarpnc
        os.system('clear')
        vars = self.get_var()
        self.plotted_variables = self.plotted_variables.union(vars)
        methods, = self.get_methods()
        self.options.logscale = bool(self.logscale.get())
        self.options.coastlines = bool(self.coastlines.get())
        self.options.countries = bool(self.countries.get())
        self.options.states = bool(self.states.get())
        self.options.counties = bool(self.counties.get())
        exec(self.pre_txt.get())
        for varkey in vars:
            figpath = eval(methods)(f = self.ncffile, varkey = varkey, options = self.options, before = self.before_txt.get(), after = self.after_txt.get())
            pncgen(getvarpnc(self.ncffile, list(vars) + ['TFLAG', 'time', 'latitude', 'longitude', 'latitude_bounds', 'longitude_bounds']), figpath + '.nc', verbose = False)
        exec(self.post_txt.get())


def StartTk(ncffile, options):
    TkApp(ncffile, options)

def pncview(f, options):
    # add a gui for plotting
    return StartTk(f, options)

def gettime(f):
    from PseudoNetCDF import PseudoNetCDFVariable
    from pylab import date2num
    from datetime import datetime
    if 'time' in f.variables:
        time = f.variables['time']
        unit = time.units
        tunit, datestr = unit.split(' since ')
        sdate = datetime.strptime(datestr.replace(' UTC', ''), '%Y-%m-%d %H:%M:%S')
        time = np.array([sdate + timedelta(**{tunit: t}) for t in time[:]])
        unit = 'time'
    elif 'TFLAG' in f.variables:
        x = f.variables['TFLAG'][:].copy()
        for axi in range(1, x.ndim - 1):
            x = x[:, 0]
        time = np.array([(datetime.strptime('%07d %06d' % (d, t), '%Y%j %H%M%S')) for d, t in x])
        unit = 'time'
    elif 'time' in f.dimensions:
        time = np.arange(len(f.dimensions['time']))
        unit = 'steps'
    elif 'TSTEP' in f.dimensions:
        time = np.arange(1, len(f.dimensions['TSTEP']) + 1)
        unit = 'steps'
    else:
        raise KeyError('No time found')
    return PseudoNetCDFVariable(None, 'time', 'f', ('time',), values = time[:], units = unit)
    
def timeseries(f, varkey, options, before = '', after = ''):
    outpath = getattr(options, 'outpath', '.')
    time = gettime(f)
    ax = pl.gca()
    var = f.variables[varkey]
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

def plot(f, varkey, options, before = '', after = ''):
    outpath = getattr(options, 'outpath', '.')
    ax = pl.gca()
    var = f.variables[varkey]
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

def presslat(f, varkey, options, before = '', after = ''):
    outpath = getattr(options, 'outpath', '.')
    ax = pl.gca()
    vert = getpresbnds(f)
    lat, latunit = getlatbnds(f)
    lat = np.append(lat.squeeze()[..., :2].mean(1), lat.squeeze()[-1, 2:].mean(0))
    var = f.variables[varkey]
    exec(before)
    ax = pl.gca()
    print varkey,
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    patches = ax.pcolor(lat, vert, var[:].squeeze(), norm = norm)
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    pl.colorbar(patches)
    ax.set_ylim(vert.max(), vert.min())
    ax.set_xlim(lat.min(), lat.max())
    fmt = 'png'
    figpath = os.path.join(outpath + '_PRESLAT_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def presslon(f, varkey, options, before = '', after = ''):
    outpath = getattr(options, 'outpath', '.')
    ax = pl.gca()
    vert = getpresbnds(f)
    lon, lonunit = getlonbnds(f)
    lon = np.append(lon.squeeze()[..., [0, 3]].mean(1), lon.squeeze()[-1, [1, 2]].mean(0))
    var = f.variables[varkey]
    exec(before)
    ax = pl.gca()
    print varkey,
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    patches = ax.pcolor(lon, vert, var[:].squeeze(), norm = norm)
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    pl.colorbar(patches)
    ax.set_ylim(vert.max(), vert.min())
    ax.set_xlim(lon.min(), lon.max())
    fmt = 'png'
    figpath = os.path.join(outpath + '_PRESLON_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath



def tileplot(f, varkey, options, before = '', after = ''):
    outpath = getattr(options, 'outpath', '.')
    ax = pl.gca()
    var = f.variables[varkey]
    exec(before)
    ax = pl.gca()
    print varkey,
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    patches = ax.pcolor(var[:].squeeze(), norm = norm)
    ax.set_xlim(0, var.squeeze().shape[1])
    ax.set_ylim(0, var.squeeze().shape[0])
    #ax.set_xlabel(X.units.strip())
    #ax.set_ylabel(Y.units.strip())
    cbar = pl.colorbar(patches)
    cbar.ax.text(.5, 1, '%.2g' % var[:].max(), horizontalalignment = 'center', verticalalignment = 'bottom')
    cbar.ax.text(.5, 0, '%.2g' % var[:].min(), horizontalalignment = 'center', verticalalignment = 'top')
    fmt = 'png'
    figpath = os.path.join(outpath + '_2D_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

def pres_from_sigma(sigma, pref, ptop, avg = False):
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres


def getsigmabnds(f):
    if hasattr(f, 'VGLVLS'):
        return f.VGLVLS[:]
    elif 'layer_edges' in f.variables:
        lay = f.variables['layer_edges']
        if lay.units.strip() in ('Pa', 'hPa'):
            sigma = (lay[:] -lay[-1]) / (lay[0] - lay[-1])
            return sigma
        else:
            warn("Unknown tranform of layer to sigma; sigma units %s" % lay.units)
            return lay
    else:
        warn("Unknown vertical coordinate")
        return np.arange(f.NLAYS)
        

def getpresmid(f):
    presb = getpresbnds(f)
    return presb[:-1] + np.diff(presb) / 2

def getsigmamid(f):
    sigmab = getsigmabnds(f)
    return sigmab[:-1] + np.diff(sigmab) / 2

def getpresbnds(f):
    if 'layer_edges' in f.variables:
        return f.variables['layer_edges'][:]
    else:
        sigma = getsigmabnds(f)
        if hasattr(f, 'VGTOP'):
            ptop = f.VGTOP
        else:
            warn("Assuming VGTOP = 100 hPa")
            ptop = 10000
        return pres_from_sigma(sigma, pref = 101325., ptop = ptop)

def getlatbnds(f):
    if 'latitude_bounds' in f.variables:
        latb = f.variables['latitude_bounds']
        unit = latb.units.strip()
        if 'nv' in latb.dimensions and latb[:].ndim == 2:
            latb = np.append(latb[:][:, 0], latb[:][-1, 1])
    elif 'latitude' in f.variables:
        latb = f.variables['latitude']
        unit = latb.units.strip()
        latdiff = np.diff(latb, axis = 0)
        if (latdiff == latdiff[0]).all():
            latb = latb - .5 * latdiff[0]
            latb = np.append(latb, latb[-1] + latdiff[0])
    elif 'ROW' in f.dimensions:
        unit = 'x (LCC km)'
        latb = np.arange(len(f.dimensions['ROW']) + 1) * getattr(f, 'YCELL', 1) / 1000.
    else:
        raise KeyError('latitude bounds not found')
    return latb, unit

def getybnds(f):
    if 'ROW' in f.dimensions:
        unit = 'x (LCC km)'
        latb = np.arange(len(f.dimensions['ROW']) + 1) * getattr(f, 'YCELL', 1)
    else:
        raise KeyError('latitude bounds not found')
    return latb, unit

def getlonbnds(f):
    if 'longitude_bounds' in f.variables:
        lonb = f.variables['longitude_bounds']
        unit = lonb.units.strip()
        if 'nv' in lonb.dimensions and lonb[:].ndim == 2:
            lonb = np.append(lonb[:][:, 0], lonb[:][-1, 1])
    elif 'longitude' in f.variables:
        lonb = f.variables['longitude']
        unit = lonb.units.strip()
        londiff = np.diff(lonb, axis = 0)
        if (londiff == londiff[0]).all():
            lonb = lonb - .5 * londiff[0]
            lonb = np.append(lonb, lonb[-1] + londiff[0])
    else:
        raise KeyError('longitude bounds not found')
    return lonb, unit

def getxbnds(f):
    if 'COL' in f.dimensions:
        unit = 'x (LCC km)'
        lonb = np.arange(len(f.dimensions['COL']) + 1) * getattr(f, 'XCELL', 1)
    else:
        raise KeyError('x bounds not found')
    return lonb, unit

def getmap(f):
    from mpl_toolkits.basemap import Basemap
    if all([hasattr(f, k) for k in 'P_GAM P_ALP P_BET XORIG YORIG XCELL YCELL'.split()]):
        llcrnrx = f.XORIG
        urcrnrx = f.XORIG + len(f.dimensions['COL']) * f.XCELL

        llcrnry = f.YORIG
        urcrnry = f.YORIG + len(f.dimensions['ROW']) * f.YCELL
        m = Basemap(projection = 'lcc', lon_0=f.P_GAM, lat_1 = f.P_ALP, lat_2 = f.P_BET, lat_0 = f.YCENT, llcrnrx = llcrnrx, llcrnry = llcrnry, urcrnry = urcrnry, urcrnrx = urcrnrx, resolution = 'i', suppress_ticks = False)
        print 'Found IO/API Mapping parameters'
    else:
        kwds = dict(suppress_ticks = False)
        try:
            lat, latunit = getlatbnds(f)
            lon, lonunit = getlonbnds(f)
            kwds['llcrnrlat'] = lat[:].min()
            kwds['urcrnrlat'] = lat[:].max()
            kwds['llcrnrlon'] = lon[:].min()
            kwds['urcrnrlon'] = lon[:].max()
        except:
            pass
        m = Basemap(**kwds)
    return m

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

def profile(f, varkey, options, before = '', after = ''):
    print varkey,
    outpath = getattr(options, 'outpath', '.')
    try:
        vert = getpresmid(f)
        vunit = 'Pa'
    except:
        vert = getsigmamid(f)
        vunit = r'\sigma'
    var = f.variables[varkey]
    dims = list(var.dimensions)
    for knownvert in ('layer', 'LAY'):
        if knownvert in dims:
            vidx = dims.index(knownvert)
            break
    else:
        raise KeyError("No known vertical coordinate; got %s" % str(dims))
    units = var.units.strip()
    vals = np.rollaxis(var[:], vidx, start = 0).reshape(vert.size, -1)
    ax = pl.gca()
    minmaxmean(ax, vals, vert)
    ax.set_xlabel(varkey + ' ('+units+')')
    ax.set_ylim(vert.max(), vert.min())
    if options.logscale: ax.set_xscale('log')
    fmt = 'png'
    figpath = os.path.join(outpath + '_map_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath

    
def mapplot(f, varkey, options, before = '', after = ''):
    outpath = getattr(options, 'outpath', '.')
    map = getmap(f)
    if map.projection == 'cyl':
        latb, latunit = getlatbnds(f)[:]
        lonb, lonunit = getlonbnds(f)[:]
    else:
        latb, latunit = getybnds(f)[:]
        lonb, lonunit = getxbnds(f)[:]
    if latb.ndim == lonb.ndim and lonb.ndim == 2:
        LON, LAT = lonb, latb
    else:
        LON, LAT = np.meshgrid(lonb, latb)

    ax = pl.gca()
    if options.logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    var = f.variables[varkey]
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
    fmt = 'png'
    figpath = os.path.join(outpath + '_map_' + varkey + '.' + fmt)
    exec(after)
    pl.savefig(figpath)
    print 'Saved fig', figpath
    return figpath
            
def main():
    from pncparse import pncparser
    ifiles, options = pncparser(has_ofile = True)
    if len(ifiles) != 1:
        raise IOError('pncview can operate on only 1 file; user requested %d' % len(ifiles))
    ifile, = ifiles
    pncview(ifile, options)

if __name__ == '__main__':
    main()
