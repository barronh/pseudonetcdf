from Tkinter import Checkbutton, Frame, Label, Scrollbar, Listbox, Button, IntVar, Tk, VERTICAL, EXTENDED, END, N, S, SINGLE, Entry, StringVar, Text, DISABLED, PhotoImage, LEFT, E, W
import os
from types import MethodType
import pylab as pl
from matplotlib.colors import Normalize, LogNorm
import numpy as np

class TkApp:
    def __init__(self, master, ncffile, options):
        self.ncffile = ncffile
        self.options = options
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
        self.pre_txt.set('pl.figure();')
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
        self.after_txt.set('ax.set_title(varkey);')
        after_label = Label(codeframe, text = 'After each figure, execute code')
        after_label.grid(row = 6, sticky = 'W')
        self.after = Entry(codeframe, width = 120, textvariable = self.after_txt)
        self.after.grid(row =7, sticky = 'E')

        self.post_txt = StringVar()
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
                
        self.methods = ['plot', 'timeseries', 'tileplot', 'mapplot']
        method_labels= ['Plot', 'Time Series', 'Tile Plot', 'Map']
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
        methods, = self.get_methods()
        self.options.logscale = bool(self.logscale.get())
        self.options.coastlines = bool(self.coastlines.get())
        self.options.countries = bool(self.countries.get())
        self.options.states = bool(self.states.get())
        self.options.counties = bool(self.counties.get())
        exec(self.pre_txt.get())
        figpath = eval(methods)(f = self.ncffile, vars = vars, options = self.options, before = self.before_txt.get(), after = self.after_txt.get())
        exec(self.post_txt.get())
        self.showfig(figpath)
    
    def showfig(self, figpath):
        import Image

        photoimage = Image.open(figpath)
        photoimage.show()


def StartTk(ncffile, options):
    root = Tk(className = 'PNCVIEW')
    
    app = TkApp(root, ncffile, options)
    
    root.mainloop()

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
        time = np.arange(len(f.dimensiosn['time']))
        unit = 'steps'
    elif 'TSTEP' in f.dimensions:
        time = np.arange(len(f.dimensiosn['TSTEP'])), 'steps'
        unit = 'steps'
    else:
        raise KeyError('No time found')
    return PseudoNetCDFVariable(None, 'time', 'f', ('time',), values = time[:], units = unit)
    
def timeseries(f, vars, options, before, after):
    outpath = getattr(options, 'outpath', '.')
    time = gettime(f)
    for varkey in vars:
        ax = pl.gca()
        exec(before)
        ax = pl.gca()
        var = f.variables[varkey]
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

def plot(f, vars, options, before, after):
    outpath = getattr(options, 'outpath', '.')
    for varkey in vars:
        ax = pl.gca()
        exec(before)
        ax = pl.gca()
        var = f.variables[varkey]
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

def tileplot(f, vars, options, before, after):
    outpath = getattr(options, 'outpath', '.')
    for varkey in vars:
        ax = pl.gca()
        exec(before)
        ax = pl.gca()
        var = f.variables[varkey]
        print varkey,
        patches = ax.pcolor(var[:].squeeze())
        #ax.set_xlabel(X.units.strip())
        #ax.set_ylabel(Y.units.strip())
        pl.colorbar(patches)
        fmt = 'png'
        figpath = os.path.join(outpath + '_2D_' + varkey + '.' + fmt)
        exec(after)
        pl.savefig(figpath)
        print 'Saved fig', figpath
    return figpath

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
    if all([hasattr(f, k) for k in 'P_GAM P_ALP P_BET XORIG YORIG'.split()]):
        llcrnrx=f.XORIG
        urcrnrx = f.XORIG + len(f.dimensions['COL']) * f.XCELL

        llcrnry=f.YORIG
        urcrnry = f.YORIG + len(f.dimensions['ROW']) * f.YCELL

        m = Basemap(projection = 'lcc', lon_0=f.P_GAM, lat_1 = f.P_ALP, lat_2 = f.P_BET, llcrnrx = llcrnrx, llcrnry = llcrnry, urcrnry = urcrnry, urcrnrx = urcrnrx, resolution = 'i', suppress_ticks = False)
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

def mapplot(f, vars, options, before, after):
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
    for varkey in vars:
        ax = pl.gca()
        if options.logscale:
            norm = LogNorm()
        else:
            norm = Normalize()
        exec(before)
        ax = pl.gca()
        var = f.variables[varkey]
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
    from pncgen import main as pncgenmain
    ofile, options = pncgenmain()
    pncview(ofile, options)

if __name__ == '__main__':
    main()