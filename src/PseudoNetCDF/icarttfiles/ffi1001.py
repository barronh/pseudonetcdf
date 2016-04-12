from __future__ import print_function
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFMaskedVariable as PseudoNetCDFVariable
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
from numpy.ma import MaskedArray, filled
from datetime import datetime, timedelta
import re
import yaml
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO
from warnings import warn

def get_lodval(v):
    try:
        return eval(v)
    except:
        return v

loddelim = re.compile('(;\s?)|,|\s')

PI_LINE = 2
ORG_LINE = 3
PLAT_LINE = 4
MISSION_LINE = 5
VOL_LINE = 6
DATE_LINE = 7
TIME_INT_LINE = 8
UNIT_LINE = 9
DATE_VAR_LINE = 10
SCALE_LINE = 11
MISSING_LINE = 12
class ffi1001(PseudoNetCDFFile):
    def __init__(self,path):
        PseudoNetCDFFile.__init__(self)
        f = open(path, 'r')
        missing = []
        units = []
        l = f.readline()
        if ',' in l:
            delim = ','
        else:
            delim = None
        split = lambda s: list(map(str.strip, s.split(delim)))

        if split(l)[-1] != '1001':
            raise TypeError("File is the wrong format.  Expected 1001; got %s" % (split(l)[-1],))
        
        n, self.fmt = split(l)
        n_user_comments = 0
        n_special_comments = 0
        self.n_header_lines = int(n)
        try:
            for li in range(self.n_header_lines-1):
                li += 2
                l = f.readline()
                LAST_VAR_DESC_LINE = 12+len(missing)
                SPECIAL_COMMENT_COUNT_LINE = LAST_VAR_DESC_LINE + 1
                LAST_SPECIAL_COMMENT_LINE = SPECIAL_COMMENT_COUNT_LINE + n_special_comments
                USER_COMMENT_COUNT_LINE = 12+len(missing)+2+n_special_comments
                if li == PI_LINE:
                    self.PI_NAME = l.strip()
                elif li == ORG_LINE:
                    self.ORGANIZATION_NAME = l.strip()
                elif li == PLAT_LINE:
                    self.SOURCE_DESCRIPTION = l.strip()
                elif li == MISSION_LINE:
                    self.MISSION_NAME = l.strip()
                elif li == VOL_LINE:
                    self.VOLUME_INFO = l.strip()
                elif li == DATE_LINE:
                    l = l.replace(',', '').split()
                    SDATE = "".join(l[:3])
                    WDATE = "".join(l[3:])
                    self.SDATE = SDATE
                    self.WDATE = WDATE
                    self._SDATE = datetime.strptime(SDATE, '%Y%m%d')
                    self._WDATE = datetime.strptime(WDATE, '%Y%m%d')
                elif li == TIME_INT_LINE:
                    self.TIME_INTERVAL = l.strip()
                elif li == UNIT_LINE:
                    units.append(l.replace('\n', '').replace('\r', '').strip())
                    self.INDEPENDENT_VARIABLE = units[-1]
                elif li == SCALE_LINE:
                    scales = [eval(i) for i in split(l)]
                    if set([float(s) for s in scales]) != set([1.]):
                        raise ValueError("Unsupported: scaling is unsupported.  data is scaled by %s" % (str(scales),))
                elif li == MISSING_LINE:
                    missing = [eval(i) for i in split(l)]
                elif li > MISSING_LINE and li <= LAST_VAR_DESC_LINE:
                    nameunit = l.replace('\n','').split(',')
                    name = nameunit[0].strip()
                    if len(nameunit) > 1:
                        units.append(nameunit[1].strip())
                    elif re.compile('(.*)\((.*)\)').match(nameunit[0]):
                        desc_groups = re.compile('(.*)\((.*)\).*').match(nameunit[0]).groups()
                        name = desc_groups[0].strip()
                        units.append(desc_groups[1].strip())
                    elif '_' in name:
                        units.append(name.split('_')[1].strip())
                    else:
                        warn('Could not find unit in string: "%s"' % l)
                        units.append(name.strip())
                elif li == SPECIAL_COMMENT_COUNT_LINE:
                    n_special_comments = int(l.replace('\n', ''))
                elif li > SPECIAL_COMMENT_COUNT_LINE and li <= LAST_SPECIAL_COMMENT_LINE:
                    pass
                elif li == USER_COMMENT_COUNT_LINE:
                    n_user_comments = int(l.replace('\n',''))
                elif li > USER_COMMENT_COUNT_LINE and li < self.n_header_lines:
                    colon_pos = l.find(':')
                    k = l[:colon_pos].strip()
                    v = l[colon_pos+1:].strip()
                    setattr(self,k,v)
                elif li == self.n_header_lines:
                    variables = l.replace(',','').split()
                    self.TFLAG = variables[0]
        except Exception as e:
            raise SyntaxError("Error parsing icartt file %s: %s" % (path, repr(e)))

        missing = missing[:1]+missing
        scales = [1.]+scales
        
        if hasattr(self,'LLOD_FLAG'):
            llod_values = loddelim.sub('\n', self.LLOD_VALUE).split()
            if len(llod_values) == 1:
                llod_values *= len(variables)
            else:
                llod_values = ['N/A']+llod_values
            
            assert len(llod_values) == len(variables)
            llod_values = [get_lodval(llod_val) for llod_val in llod_values]
            
            llod_flags = len(llod_values)*[self.LLOD_FLAG]
            llod_flags = [get_lodval(llod_flag) for llod_flag in llod_flags]
        
        if hasattr(self,'ULOD_FLAG'):
            ulod_values = loddelim.sub('\n', self.ULOD_VALUE).split()
            if len(ulod_values) == 1:
                ulod_values *= len(variables)
            else:
                ulod_values = ['N/A']+ulod_values

            assert len(ulod_values) == len(variables)
            ulod_values = [get_lodval(ulod_val) for ulod_val in ulod_values]
            
            ulod_flags = len(ulod_values)*[self.ULOD_FLAG]
            ulod_flags = [get_lodval(ulod_flag) for ulod_flag in ulod_flags]
        
        data = f.read()
        datalines = data.split('\n')
        ndatalines = len(datalines)
        while datalines[-1] in ('', ' ', '\r'):
            ndatalines -=1
            datalines.pop(-1)
        data = genfromtxt(StringIO(u'\n'.join(datalines)), delimiter = delim, dtype = 'd')
        data = data.reshape(ndatalines,len(variables))
        data = data.swapaxes(0,1)
        self.createDimension('POINTS', ndatalines)
        for var, scale, miss, unit, dat, llod_flag, llod_val, ulod_flag, ulod_val in zip(variables, scales, missing, units, data, llod_flags, llod_values, ulod_flags, ulod_values):
            vals = MaskedArray(dat, mask = dat == miss, fill_value = miss)
            tmpvar = self.variables[var] = PseudoNetCDFVariable(self, var, 'd', ('POINTS',), values = vals)
            tmpvar.units = unit
            tmpvar.standard_name = var
            tmpvar.missing_value = miss
            tmpvar.fill_value = miss
            tmpvar.scale = scale

            if hasattr(self,'LLOD_FLAG'):
                tmpvar.llod_flag = llod_flag
                tmpvar.llod_value = llod_val

            if hasattr(self,'ULOD_FLAG'):
                tmpvar.ulod_flag = ulod_flag
                tmpvar.ulod_value = ulod_val

        
        self._date_objs = self._SDATE + vectorize(lambda s: timedelta(seconds = int(s), microseconds = (s - int(s)) * 1.E6 ))(self.variables[self.TFLAG]).view(type = ndarray)

def ncf2ffi1001(f, outpath, mode = 'w'):
    outfile = open(outpath, mode)
    header_keys = "PI_CONTACT_INFO PLATFORM LOCATION ASSOCIATED_DATA INSTRUMENT_INFO DATA_INFO UNCERTAINTY ULOD_FLAG ULOD_VALUE LLOD_FLAG LLOD_VALUE DM_CONTACT_INFO PROJECT_INFO STIPULATIONS_ON_USE OTHER_COMMENTS REVISION".split()
    print('%d, %d' % (len(f.ncattrs()) + len(f.variables), 1001), file = outfile)
    print(getattr(f, 'PI_NAME', 'Unknown'), file = outfile)
    print(getattr(f, 'ORGANIZATION_NAME', 'Unknown'), file = outfile)
    print(getattr(f, 'SOURCE_DESCRIPTION', 'Unknown'), file = outfile)
    print(getattr(f, 'VOLUME_INFO', 'Unknown'), file = outfile)
    print(f.SDATE, f.WDATE, file = outfile)
    print(f.TIME_INTERVAL, file = outfile)
    print(f.INDEPENDENT_VARIABLE, file = outfile)
    print('%d' % len(f.variables), file = outfile)
    for key, var in f.variables.items():
        print('%s, %s' % (key, getattr(var, 'units', 'unknown')), file = outfile)
    
    print(len(f.ncattrs()), file = outfile)
    for key in f.ncattrs():
        print('%s: %s' % (key, getattr(f, key, '')), file = outfile)
    
    vals = [filled(f.variables[f.INDEPENDENT_VARIABLE][:]).ravel()]
    keys = [f.INDEPENDENT_VARIABLE]
    for key, var in f.variables.items():
        if key == f.INDEPENDENT_VARIABLE: continue
        keys.append(key)
        vals.append(filled(var[:]).ravel())
        
    print(', '.join(keys), file = outfile)
    for row in array(vals).T:
        row.tofile(outfile, format = '%.6e', sep = ', ')
        print('', file = outfile)

from PseudoNetCDF._getwriter import registerwriter
registerwriter('ffi1001', ncf2ffi1001)    
