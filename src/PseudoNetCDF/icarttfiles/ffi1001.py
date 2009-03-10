from ..sci_var import PseudoNetCDFFile, PseudoNetCDFVariable
from numpy import fromstring, vectorize, ndarray, array
from numpy.ma import MaskedArray
from datetime import datetime, timedelta
import re
from warnings import warn

class ffi1001(PseudoNetCDFFile):
    def __init__(self,path):
        self.dimensions = {}
        self.variables = {}
        f = file(path, 'r')
        missing = []
        units = []
        l = f.readline()
        if l.split()[-1] != '1001':
            raise TypeError, "File is the wrong format.  Expected 1001; got %s" % (l.split()[-1],)
        
        n, self.fmt = l.split()
        self.n_header_lines = int(n)
        for li in range(self.n_header_lines-1):
            li += 2
            l = f.readline()
            if li == 7:
                l = l.replace(',', '').split()
                SDATE = " ".join(l[:3])
                WDATE = " ".join(l[3:])
                self.SDATE = datetime.strptime(SDATE, '%Y %m %d')
                self.WDATE = datetime.strptime(WDATE, '%Y %m %d')
            elif li == 9:
                units.append(l.replace('\n', ''))
            elif li == 11:
                scales = l.split()
                if set([float(s) for s in scales]) != set([1.]):
                    raise ValueError, "Unsupported: scaling is unsupported.  data is scaled by %s" % (str(scales),)
            elif li == 12:
                missing = l.split()
            elif li > 12 and li <= 12+len(missing):
                nameunit = l.replace('\n','').split(',')
                name = nameunit[0].strip()
                if len(nameunit) > 1:
                    units.append(nameunit[1])
                elif re.compile('(.*)\((.*)\)').match(nameunit[0]):
                    desc_groups = re.compile('(.*)\((.*)\).*').match(nameunit[0]).groups()
                    name = desc_groups[0].strip()
                    units.append(desc_groups[1].strip())
                elif '_' in name:
                    units.append(name.split('_')[1])
                else:
                    warn('Could not find unit in string: "%s"' % l)
                    units.append(name)

            elif li == self.n_header_lines:
                variables = l.replace(',','').split()
                self.TFLAG = variables[0]

        missing = missing[:1]+missing
        scales = [1.]+scales
        data = f.read()
        datalines = data.split('\n')
        ndatalines = len(datalines)
        while datalines[-1] in ('', ' ', '\r'):
            ndatalines -=1
            datalines.pop(-1)
        data = fromstring(data,dtype ='d', sep = ' ')
        data = data.reshape(ndatalines,len(variables))
        data = data.swapaxes(0,1)
        self.createDimension('POINTS', ndatalines)
        for var, scale, miss, unit, dat in zip(variables, scales, missing, units, data):
            vals = MaskedArray(dat, mask = dat == miss)
            tmpvar = self.variables[var] = PseudoNetCDFVariable(self, var, 'f', ('POINTS',), values = vals)
            tmpvar.units = unit

            tmpvar.fillvalue = miss
            tmpvar.scale = scale
            tmpvar[:] = dat
        
        
        self.date_objs = self.SDATE + vectorize(lambda s: timedelta(seconds = int(s), microseconds = (s - int(s)) * 1.E6 ))(self.variables[self.TFLAG]).view(type = ndarray)
        self.createDimension('YYYYMMDDTHHMMSS.microS', 22)
        var = self.createVariable('TFLAG', 'c', ('POINTS', 'YYYYMMDDTHHMMSS.microS'))
        var[:] = array([d.strftime('%Y%m%dT%H%M%S.%f') for d in self.date_objs], dtype = '|S22').view('|S1').reshape(self.date_objs.shape[0], self.dimensions['YYYYMMDDTHHMMSS.microS'])
        var.units = 'YYYYMMDDTHHMMSS.microS'
        var.fillvalue = ''
        var.scale = 1.
