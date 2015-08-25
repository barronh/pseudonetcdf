__all__ = ['cdtoms']

from PseudoNetCDF import PseudoNetCDFFile
from re import compile
from numpy import array, arange
from datetime import datetime
dayre = compile(' Day:\s+(?P<jday>\d+) (?P<daystring>.{12})\s+EP/TOMS CORRECTED OZONE GEN:\d+\.\d+\sV\d ALECT:\s+\d+:\d+ [AP]M ')
lonre = compile(' Longitudes:\s+(?P<lonbins>\d+)\s+bins centered on\s+(?P<blon>\d+\.\d+\s+[WE]) to\s+(?P<elon>\d+\.\d+\s+[WE])\s+\((?P<lonbinsize>\d+.\d+) degree steps\)  ')
latre = compile(' Latitudes\s:\s+(?P<latbins>\d+)\s+bins centered on\s+(?P<blat>\d+\.\d+\s+[SN]) to\s+(?P<elat>\d+\.\d+\s+[SN])\s+\((?P<latbinsize>\d+.\d+) degree steps\)  ')


def cdtoms(path):
    outfile = PseudoNetCDFFile()
    inlines = open(path, 'rU').readlines()
    dayline = dayre.match(inlines[0]).groupdict()
    date = datetime.strptime(dayline['daystring'], '%b %d, %Y')
    lonline = lonre.match(inlines[1]).groupdict()
    latline = latre.match(inlines[2]).groupdict()
    for propdict in [dayline, lonline, latline]:
        for k,v in propdict.iteritems():
            try:
                v = eval(v)
            except:
                pass
            setattr(outfile, k, v)
    blat, bsn = outfile.blat.split()
    elat, esn = outfile.elat.split()
    blat = {'N': 1, 'S': -1}[bsn] * eval(blat)
    elat = {'N': 1, 'S': -1}[esn] * eval(elat)

    blon, bwe = outfile.blon.split()
    elon, ewe = outfile.elon.split()
    blon = {'E': 1, 'W': -1}[bwe] * eval(blon)
    elon = {'E': 1, 'W': -1}[ewe] * eval(elon)


    outfile.createDimension('LAT', outfile.latbins)
    outfile.createDimension('LON', outfile.lonbins)
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
    var = outfile.createVariable('ozone', 'f', ('LAT', 'LON'))
    var.units = 'matm-cm'
    var.long_name = var.var_desc = 'ozone'.ljust(16)
    var[:] = array([[eval(s[n*i:n*i+n]) for i in range(len(s)/n)] for s in datalines if s.strip() != ''], 'f')

    var = outfile.createVariable('lat', 'f', ('LAT',))
    var.units = 'degrees N'
    var[:] = arange(blat, elat+outfile.latbinsize, outfile.latbinsize)

    var = outfile.createVariable('lon', 'f', ('LON',))
    var.units = 'degrees E'
    var[:] = arange(blon, elon+outfile.lonbinsize, outfile.lonbinsize)


    return outfile
    
if __name__ == '__main__':
    from PseudoNetCDF.pncdump import pncdump
    pfile = cdtoms('test.txt')
    pncdump(pfile)
    print pfile.variables['ozone'][2,3]
    