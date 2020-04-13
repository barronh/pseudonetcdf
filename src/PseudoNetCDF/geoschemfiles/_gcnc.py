__all__ = ['gcnc_base', 'gcnc']
from PseudoNetCDF.core._files import netcdf
from PseudoNetCDF import PseudoNetCDFFile
from PseudoNetCDF.pncwarn import warn
import numpy as np


class gcnc_base(PseudoNetCDFFile):
    def _newlike(self):
        if isinstance(self, PseudoNetCDFFile):
            outt = gcnc_base
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        return outf

    def ll2xy(self, lon, lat, *args, **kwds):
        return lon, lat

    def ll2ij(self, lon, lat, bounds='error', clean='none'):
        """
        Converts lon/lat to 0-based indicies (0,M), (0,N)

        Parameters
        ----------
        lon : scalar or iterable of longitudes in decimal degrees
        lat : scalar or iterable of latitudes in decimal degrees
        bounds : ignore, error, warn if i,j are out of domain

        Returns
        -------
        i, j : indices (0-based) for variables
        """
        lon = np.asarray(lon)
        lat = np.asarray(lat)

        lonc = self.variables['lon']
        latc = self.variables['lat']

        dlon = np.median(np.diff(lonc)) / 2.
        dlat = np.median(np.diff(latc)) / 2.

        lonb = np.array([lonc - dlon, lonc + dlon]).T
        latb = np.array([latc - dlat, latc + dlat]).T

        spi = latb[:, 0] < -90
        npi = latb[:, 0] > 90

        latb[spi, 0] = latc[spi] - dlon / 2.
        latb[spi, 1] = latc[spi] + dlon / 2.
        latb[npi, 0] = latc[npi] - dlon / 2.
        latb[npi, 1] = latc[npi] + dlon / 2.

        easter = lon[:, None] >= lonb[:, 0]
        wester = lon[:, None] < lonb[:, 1]
        loni, gi = np.where(easter & wester)
        i = np.ma.masked_all(lon.shape, dtype='i')
        i[loni] = gi
        norther = lat[:, None] >= latb[:, 0]
        souther = lat[:, None] < latb[:, 1]
        latj, gj = np.where(norther & souther)
        j = np.ma.masked_all(lat.shape, dtype='i')
        j[latj] = gj

        nx = lonb.shape[0]
        ny = latb.shape[0]
        if bounds == 'ignore':
            pass
        else:
            missi, = np.where(i.mask)
            missj, = np.where(j.mask)
            lowi, = np.where(~easter[missi].any(1))
            highi, = np.where(~wester[missi].any(1))
            lowj, = np.where(~norther[missj].any(1))
            highj, = np.where(~souther[missj].any(1))
            outb = (lowi | lowj | highi | highj)
            nout = outb.sum()
            if nout > 0:
                message = '{} Points out of bounds; {}'.format(
                    nout, np.where(outb))
                if bounds == 'error':
                    raise ValueError(message)
                else:
                    warn(message)

        if clean == 'clip':
            i[lowi] = 0
            i[highi] = nx - 1
            j[lowj] = 0
            j[highj] = ny - 1
        else:
            # Mask or nothing should both create symmetric masks
            mask = (
                np.ma.getmaskarray(i) |
                np.ma.getmaskarray(j)
            )
            i = np.ma.masked_where(mask, i)
            j = np.ma.masked_where(mask, j)

        return i, j

    def interpSigma(self, vglvls, vgtop=None, interptype='linear',
                    extrapolate=False, fill_value='extrapolate',
                    layerdims=None,
                    approach='eta', verbose=0):
        """
        Parameters
        ----------
        self : the file to interpolate from must have VGLVLS
        vglvls : the new vglvls (edges)
        vgtop : Converting to new vgtop (Pascals)
        interptype : 'linear' or 'conserve'
             linear : uses a linear interpolation
             conserve : uses a mass conserving interpolation
        extrapolate : allow extrapolation beyond bounds with linear,
                      default False
        fill_value : set fill value (e.g, nan) to prevent extrapolation or edge
                     continuation
        layerdims : specify layer dimension, None will apply to all dimensions
                    named lev*
        approach :
             eta : use simple eta coordinates to calculate sigma and
                   interpolate
             pressure : requires surface pressure
        verbose : 0-inf show more

        Returns
        -------
        outf - ioapi_base PseudoNetCDFFile with al variables interpolated

        Notes
        -----
        When extrapolate is false, the edge values are used for points beyond
        the inputs.
        """
        etai = (
            self.variables['P0'][...] * self.variables['hybi'][...] +
            self.variables['hyai'][...]
        ) * 100.
        # grab input sigma coordinates
        myvglvls = (etai - vgtop) / (etai[0] - vgtop)

        # Use midpoint for sigmas in inputs and outputs
        zs = (myvglvls[:-1] + myvglvls[1:]) / 2.
        nzs = (vglvls[:-1] + vglvls[1:]) / 2.

        if interptype == 'linear':
            from ..coordutil import getinterpweights
            weights = getinterpweights(
                zs, nzs, kind=interptype, fill_value=fill_value,
                extrapolate=extrapolate)
            # Create a function for interpolation

            def interpsigma(data):
                if data.ndim == 1:
                    newdata = (weights[:data.shape[0]] * data[:, None]).sum(0)
                else:
                    newdata = (weights[None, :data.shape[0], :, None, None] *
                               data[:, :, None]).sum(1)
                return newdata
        else:
            raise ValueError(
                'interptype only implemented for "linear"; got ' + interptype
            )

        # Apply function on LAY
        if layerdims is None:
            layerdims = sorted([
                dk for dk in self.dimensions
                if (
                    dk.startswith('lev') and
                    not dk.endswith('_bounds') and
                    not dk == 'lev1')
            ])
            if verbose > 0:
                print(layerdims)
        dimfunc = dict([(layerkey, interpsigma) for layerkey in layerdims])
        outf = self.applyAlongDimensions(**dimfunc)

        return outf

    def stack(self, other, dimkey):
        from collections import Iterable
        outf = self._copywith(props=True, dimensions=False)
        if isinstance(other, Iterable):
            pfiles = other
        else:
            pfiles = [other]

        outf = PseudoNetCDFFile.stack(self, other, dimkey)
        tvar = outf.variables['time']
        tunits = tvar.units.strip()
        tres, rdatestr = tunits.split(' since ')

        if dimkey == 'time':
            times = self.getTimes()
            for pfile in pfiles:
                times = np.concatenate([times, pfile.getTimes()], axis=0)

        tformat = '%Y-%m-%d %H:%M:%S%z'
        rdate = times[0]
        rdatestr = rdate.strftime(tformat)
        tvar.units = 'hours since ' + rdatestr
        rvals = np.array([
            (t - rdate).total_seconds() / 3600. for t in times
        ])
        tvar[:] = rvals
        return outf


class gcnc(gcnc_base, netcdf):
    def __init__(self, *args, **kwds):
        netcdf.__init__(self, *args, **kwds)
        self.setCoords('hyam hyai hybm hybi time lat lon lev ilev P0'.split())
