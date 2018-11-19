from __future__ import print_function
__all__ = ['cdtoms']


from PseudoNetCDF import PseudoNetCDFFile
from re import compile
from numpy import array, arange


# from datetime import datetime
dayre = compile(
    r' Day:\s+(?P<jday>\d+) (?P<daystring>.{12})\s+EP/' +
    r'TOMS CORRECTED OZONE GEN:\d+\.\d+\sV\d ALECT:\s+\d+' +
    r':\d+ [AP]M '
)
lonre = compile(
    r' Longitudes:\s+(?P<lonbins>\d+)\s+bins centered ' +
    r'on\s+(?P<blon>\d+\.\d+\s+[WE]) to\s+(?P<elon>\d+\.' +
    r'\d+\s+[WE])\s+\((?P<lonbinsize>\d+.\d+) degree steps\)  '
)
latre = compile(
    r' Latitudes\s:\s+(?P<latbins>\d+)\s+bins centered on' +
    r'\s+(?P<blat>\d+\.\d+\s+[SN]) to\s+(?P<elat>\d+\.\d+\s+[SN])' +
    r'\s+\((?P<latbinsize>\d+.\d+) degree steps\)  '
)


class cdtoms(PseudoNetCDFFile):
    def __init__(self, path):
        inlines = open(path, 'rU').readlines()
        dayline = dayre.match(inlines[0]).groupdict()
        # date = datetime.strptime(dayline['daystring'], '%b %d, %Y')
        lonline = lonre.match(inlines[1]).groupdict()
        latline = latre.match(inlines[2]).groupdict()
        for propdict in [dayline, lonline, latline]:
            for k, v in propdict.items():
                try:
                    v = eval(v)
                except Exception:
                    pass
                setattr(self, k, v)
        blat, bsn = self.blat.split()
        elat, esn = self.elat.split()
        blat = {'N': 1, 'S': -1}[bsn] * eval(blat)
        elat = {'N': 1, 'S': -1}[esn] * eval(elat)

        blon, bwe = self.blon.split()
        elon, ewe = self.elon.split()
        blon = {'E': 1, 'W': -1}[bwe] * eval(blon)
        elon = {'E': 1, 'W': -1}[ewe] * eval(elon)

        self.createDimension('LAT', self.latbins)
        self.createDimension('LON', self.lonbins)
        datalines = inlines[3:]
        lats = []
        for i, line in enumerate(datalines):
            if 'lat' not in line:
                datalines[i] = line[1:-1]
            else:
                data, lat = line.split('lat =')
                datalines[i] = data[1:-1].rstrip() + '\n'
                lats.append(lat.strip())

        n = 3

        datalines = ''.join(datalines).split('\n')
        var = self.createVariable('ozone', 'f', ('LAT', 'LON'))
        var.units = 'matm-cm'
        var.long_name = var.var_desc = 'ozone'.ljust(16)
        var[:] = array([[eval(s[n * i:n * i + n]) for i in range(int(len(s) / n))]
                        for s in datalines if s.strip() != ''], 'f')

        var = self.createVariable('lat', 'f', ('LAT',))
        var.units = 'degrees N'
        var[:] = arange(blat, elat + self.latbinsize, self.latbinsize)

        var = self.createVariable('lon', 'f', ('LON',))
        var.units = 'degrees E'
        var[:] = arange(blon, elon + self.lonbinsize, self.lonbinsize)


if __name__ == '__main__':
    from PseudoNetCDF.pncdump import pncdump
    pfile = cdtoms('test.txt')
    pncdump(pfile)
    print(pfile.variables['ozone'][2, 3])
