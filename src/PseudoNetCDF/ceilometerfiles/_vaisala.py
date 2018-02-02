from PseudoNetCDF import PseudoNetCDFFile
import numpy as np
import matplotlib.pyplot as plt

class ceilometerl2(PseudoNetCDFFile):
    """
    ceilometerhis is designed to read the "*.his"' files from
    the Vaisala ceilometers (CL31 and CL51), which is a level
    2 raw output
    """
    def __init__(self, path):
        """
        Parameters
        ----------
        path : path to output file
        """
        try:
            import pandas as pd
        except:
            raise ImportError('ceilometerl2 requires pandas; install pandas (e.g., pip install pandas)')
        self._path = path
        self._data = pd.read_csv(path, skiprows = 1, parse_dates = ['CREATEDATE']).rename(columns = lambda x: x.strip())
        self._data['BS_PROFILE'] = self._data['BS_PROFILE'].str.strip()
        nchar = len(self._data['BS_PROFILE'][0])
        nlays = nchar / 5
        data = np.array(''.join(self._data['BS_PROFILE']), dtype = 'c').reshape(-1, nlays, 5)
        BS_PROFILE = (np.vectorize(int)(data, 16) * 16**np.arange(5)[::-1]).sum(2)
        #BS_PROFILE = np.vectorize(lambda x: eval(b'0x' + x))(data.view('S5')[:,:,0])
        self.createDimension('time', data.shape[0])
        self.createDimension('altitude_tops', nlays)
        self.createDimension('altitude_edges', nlays + 1)
        var = self.createVariable('time', 'd', ('time',))
        var[:] =self._data['UNIXTIME']
        var.description = 'time'
        var.units = 'seconds since 1970-01-01 00:00:00+0000'
        var = self.createVariable('PERIOD', 'i', ('time',))
        var[:] = self._data['PERIOD']
        var.units = 'seconds'
        altitude_tops = self.createVariable('altitude_tops', 'd', ('altitude_tops',))
        altitude_tops.description = "Profile bin altitude, meters above ground level" ;
        altitude_tops.units = "meters" ;
        altitude_tops.valid_range = np.array([0, 4500])
        altitude_tops[:] = np.arange(10, 4510, 10)
        altitude_edges = self.createVariable('altitude_edges', 'd', ('altitude_edges',))
        altitude_edges.units = "meters" ;
        altitude_edges.valid_range = np.array([0, 4500])
        altitude_edges.description = "Profile bin altitude, meters above ground level" ;
        altitude_edges[:] = np.arange(0, 4510, 10)
        blview_backscatter_profile = self.createVariable('blview_backscatter_profile', 'd', ('time','altitude_tops'), fill_value = -999)
        blview_backscatter_profile.description = "Backscatter profile at 910 nm not more than 4,500 m above ground level" ;
        blview_backscatter_profile.units = "10^-9/m/sr" ;
        blview_backscatter_profile.valid_range = np.array([10, 540000])
        blview_backscatter_profile[:]  = np.ma.masked_greater(np.ma.masked_less(BS_PROFILE, 10), 540000)
    
    def plot(self, fig_kw = {'figsize': (12, 6)}, ax_kw = {}, cax_kw = {}, pcolor_kw = {}):
        """
        Parameters
        ----------
        fig_kw : keywords for figure creation, default {'figsize': (12, 6)}
        ax_kw  : keywords for axes creation of the primary axis, default {}
        cax_kw  : keywords for axes creation of the colorbar axis, default {}
        pcolor_kw : keywords for axes creation of the pseudo-color mesh, default {}.
                    the default value for norm is added as matplotlib.colors.LogNorm()
                    if not provided
        
        Retruns
        -------
        fig : figure object produced
        
        Notes
        -----
        - fig_kw : see matplotlilb.pyplot.figure for options
        - ax_kw : and cax_kw see matptlotlib.pyplot.axes for options
                  except that ax_kw and cax_kw both can take an optional rect
                  keyword the rect will be used as the positional argument to
                  figure.add_axes
          - the default primary axis rect is [.075, .15, .8, .8]
          - the default colorbar axis rect starts 1% figure width to the
            right of the primary axis and takes 1/3 of the remaining space
        - pcolor_kw : see matplotlib.pyplot.pcolormesh for options
        """
        timedelta = plt.matplotlib.dates.datetime.timedelta
        fig = plt.figure(**fig_kw)
        axrect = ax_kw.pop('rect', [.075, .15, .8, .8])
        ax = fig.add_axes(axrect, **ax_kw)
        # The default color bar axis
        caxrect = axrect.copy()
        lax = axrect[0]
        wax = axrect[2]
        margin = wax * .01
        clax = lax + wax + margin
        cwax = (1 - clax) * .3
        caxrect[0] = clax
        caxrect[2] = cwax
        caxrect = cax_kw.pop('rect', caxrect)
        cax = fig.add_axes(caxrect, **cax_kw)
        # Get times as date/time object
        times = self.getTimes()
        # Get duration of sample
        periods = self.variables['PERIOD']

        # add ending bound
        time_edges = np.append(times, times[-1] + timedelta(seconds = float(periods[-1])))

        # Get altitude edges
        altedges = self.variables['altitude_edges']
        pr = self.variables['blview_backscatter_profile']
        pcolor_kw.setdefault('norm', plt.matplotlib.colors.LogNorm())
        patches = ax.pcolormesh(time_edges[:], altedges[:], pr[:].T, **pcolor_kw)
        plt.colorbar(patches, cax = cax, label = pr.units.strip())
        return fig

if __name__ == '__main__':
    from PseudoNetCDF import pncopen
    f = pncopen('/work/ROMO/users/kpc/flintfx/convceilometer/CEILOMETER_1_LEVEL_2_20.his', format = 'ceilometerl2')
