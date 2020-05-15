from __future__ import print_function
import sys
import unittest
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from PseudoNetCDF.sci_var import PseudoNetCDFMaskedVariable
from PseudoNetCDF._getwriter import registerwriter
from numpy import vectorize, ndarray, array, genfromtxt
from numpy.ma import MaskedArray, filled
import numpy as np
from datetime import datetime, timedelta
import re
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO
from warnings import warn

if (sys.version_info > (3, 0)):
    openf = open
else:
    def openf(path, mode, encoding):
        return open(path, mode)

PseudoNetCDFVariable = PseudoNetCDFMaskedVariable


def get_lodval(v):
    try:
        return eval(v)
    except Exception:
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
    """
Overview:
  ffi1001 is a reader object for the NASA AMES format also
  known as the ICARTT file format. The format is defined in detail
  at https://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm
  Standard of practice is to write files in UTF-8 encoding.
  However, it is not uncommon to receive files with special
  characters. To specify an encoding use the encoding keyword.

Example:
  datafile = ffi1001(path, encoding = 'latin1')
    """
    @classmethod
    def isMine(cls, path):
        try:
            ffifmt = openf(
                path, 'rU', encoding='utf-8').readline().strip()[-4:]
            if hasattr(ffifmt, 'decode'):
                ffifmt = ffifmt.decode()
            testfmt = u'1001'
            return ffifmt == testfmt
        except Exception:
            return False

    def __init__(self, path, keysubs={'/': '_'}, encoding='utf-8',
                 default_llod_flag=-8888, default_llod_value='N/A',
                 default_ulod_flag=-7777, default_ulod_value='N/A'):
        """
Arguments:
   self - implied input (not supplied in call)
   path - path to file
   keysubs - dictionary of characters to remove from variable keys and
             their replacements
   encoding - file encoding (utf-8, latin1, cp1252, etc.)
   default_llod_flag - flag value for lower limit of detections if not
                       specified
   default_llod_value - default value to use for replacement of llod_flag
   default_ulod_flag - flag value for upper limit of detections if not
                       specified
   default_ulod_value - default value to use for replacement of ulod_flag
Returns:
   out - PseudoNetCDFFile interface to data in file.
        """
        lastattr = None
        PseudoNetCDFFile.__init__(self)
        f = openf(path, 'rU', encoding=encoding)
        missing = []
        units = []
        line = f.readline()
        if ',' in line:
            delim = ','
        else:
            delim = None

        def split(s):
            return [s_.strip() for s_ in s.split(delim)]

        if split(line)[-1] != '1001':
            raise TypeError(
                "File is the wrong format.  " +
                "Expected 1001; got %s" % (split(line)[-1],))

        n, self.fmt = split(line)
        # n_user_comments = 0
        n_special_comments = 0
        self.n_header_lines = int(n)
        try:
            for li in range(self.n_header_lines - 1):
                li += 2
                line = f.readline()
                LAST_VAR_DESC_LINE = 12 + len(missing)
                SPECIAL_COMMENT_COUNT_LINE = LAST_VAR_DESC_LINE + 1
                LAST_SPECIAL_COMMENT_LINE = (SPECIAL_COMMENT_COUNT_LINE +
                                             n_special_comments)
                USER_COMMENT_COUNT_LINE = (12 + len(missing) + 2 +
                                           n_special_comments)
                if li == PI_LINE:
                    self.PI_NAME = line.strip()
                elif li == ORG_LINE:
                    self.ORGANIZATION_NAME = line.strip()
                elif li == PLAT_LINE:
                    self.SOURCE_DESCRIPTION = line.strip()
                elif li == MISSION_LINE:
                    self.MISSION_NAME = line.strip()
                elif li == VOL_LINE:
                    self.VOLUME_INFO = ', '.join(split(line))
                elif li == DATE_LINE:
                    line = line.replace(',', ' ').replace(
                        '-', ' ').replace('  ', ' ').split()
                    SDATE = ", ".join(line[:3])
                    WDATE = ", ".join(line[3:])
                    self.SDATE = SDATE
                    self.WDATE = WDATE
                    self._SDATE = datetime.strptime(SDATE, '%Y, %m, %d')
                    self._WDATE = datetime.strptime(WDATE, '%Y, %m, %d')
                elif li == TIME_INT_LINE:
                    self.TIME_INTERVAL = line.strip()
                elif li == UNIT_LINE:
                    unitstr = line.replace('\n', '').replace('\r', '').strip()
                    self.INDEPENDENT_VARIABLE_DEFINITION = unitstr
                    unitstr = [s.strip() for s in unitstr.split(',')]
                    self.INDEPENDENT_VARIABLE = unitstr[0]
                    nstr = len(unitstr)
                    if nstr==1:
                        units.append(unitstr[0])
                        self.INDEPENDENT_VARIABLE_UNITS = unitstr[0]
                    if nstr>=2:
                        units.append(unitstr[1])
                        self.INDEPENDENT_VARIABLE_UNITS = unitstr[1]
                elif li == SCALE_LINE:
                    scales = [eval(i) for i in split(line)]
                elif li == MISSING_LINE:
                    missing = [eval(i) for i in split(line)]
                elif li > MISSING_LINE and li <= LAST_VAR_DESC_LINE:
                    nameunit = line.replace('\n', '').split(',')
                    name = nameunit[0].strip()
                    if len(nameunit) > 1:
                        units.append(nameunit[1].strip())
                    elif re.compile('(.*)\((.*)\)').match(nameunit[0]):
                        desc_groups = re.compile(
                            '(.*)\((.*)\).*').match(nameunit[0]).groups()
                        name = desc_groups[0].strip()
                        units.append(desc_groups[1].strip())
                    elif '_' in name:
                        units.append(name.split('_')[1].strip())
                    else:
                        warn('Could not find unit in string: "%s"' % line)
                        units.append(name.strip())
                elif li == SPECIAL_COMMENT_COUNT_LINE:
                    n_special_comments = int(line.replace('\n', ''))
                elif (li > SPECIAL_COMMENT_COUNT_LINE and
                      li <= LAST_SPECIAL_COMMENT_LINE):
                    colon_pos = line.find(':')
                    if (li==(SPECIAL_COMMENT_COUNT_LINE+1) and colon_pos==-1):
                        k = 'SPECIAL_COMMENTS'
                        v = line.strip()
                    elif (line[:1] == ' ' or colon_pos==-1):
                        # Append to prior attribute line
                        k = lastattr
                        v = getattr(self, k, '') + line.rstrip()
                    else:
                        k = line[:colon_pos].strip().replace('/','_')
                        v = line[colon_pos + 1:].strip()
                    setattr(self, k, v)
                    lastattr = k
                elif li == USER_COMMENT_COUNT_LINE:
                    lastattr = None
                    # n_user_comments = int(line.replace('\n', ''))
                elif (li > USER_COMMENT_COUNT_LINE and
                      li < self.n_header_lines):
                    colon_pos = line.find(':')
                    if line[:1] == ' ':
                        k = lastattr
                        v = getattr(self, k, '') + line
                    else:
                        k = line[:colon_pos].strip()
                        v = line[colon_pos + 1:].strip()
                    setattr(self, k, v)
                    lastattr = k
                elif li == self.n_header_lines:
                    varstr = line.replace(',', ' ').replace('  ', ' ')
                    variables = varstr.split()
                    for oc, nc in keysubs.items():
                        variables = [vn.replace(oc, nc) for vn in variables]
                    self.TFLAG = variables[0]
        except Exception as e:
            raise SyntaxError(
                "Error parsing icartt file %s: %s" % (path, repr(e)))

        missing = missing[:1] + missing
        scales = [1.] + scales

        nvars = len(variables)
        if hasattr(self, 'LLOD_FLAG'):
            llod_values = loddelim.sub('\n', self.LLOD_VALUE).split()
            if len(llod_values) == 1:
                llod_values *= nvars
            elif len(llod_values) == (nvars - 1):
                llod_values = ['N/A'] + llod_values
            else:
                warn(
                    'Expected 1 or %d LLOD_VALUE(s); got %d' %
                    (nvars - 1, len(llod_values))
                )
                llod_values = ['N/A'] * nvars

            assert len(llod_values) == len(variables)
            llod_values = [get_lodval(llod_val) for llod_val in llod_values]

            llod_flags = len(llod_values) * [self.LLOD_FLAG]
            llod_flags = [get_lodval(llod_flag) for llod_flag in llod_flags]
        else:
            llod_flags = [default_llod_flag] * len(scales)
            llod_values = [default_llod_value] * len(scales)

        if hasattr(self, 'ULOD_FLAG'):
            ulod_values = loddelim.sub('\n', self.ULOD_VALUE).split()
            if len(ulod_values) == 1:
                ulod_values *= nvars
            elif len(ulod_values) == (nvars - 1):
                ulod_values = ['N/A'] + ulod_values
            else:
                warn(
                    'Expected 1 or %d ULOD_VALUE(s); got %d' %
                    (nvars - 1, len(ulod_values))
                )
                ulod_values = ['N/A'] * nvars

            assert len(ulod_values) == len(variables)
            ulod_values = [get_lodval(ulod_val) for ulod_val in ulod_values]

            ulod_flags = len(ulod_values) * [self.ULOD_FLAG]
            ulod_flags = [get_lodval(ulod_flag) for ulod_flag in ulod_flags]
        else:
            ulod_flags = [default_ulod_flag] * len(scales)
            ulod_values = [default_ulod_value] * len(scales)

        data = f.read()
        datalines = data.split('\n')
        ndatalines = len(datalines)
        while datalines[-1] in ('', ' ', '\r'):
            ndatalines -= 1
            datalines.pop(-1)

        data = genfromtxt(StringIO('\n'.join(datalines).encode()),
                          delimiter=delim, dtype='d')
        data = data.reshape(ndatalines, len(variables))
        data = data.swapaxes(0, 1)
        self.createDimension('POINTS', ndatalines)
        for vi, var in enumerate(variables):
            scale = scales[vi]
            miss = missing[vi]
            unit = units[vi]
            dat = data[vi]
            llod_flag = llod_flags[vi]
            llod_val = llod_values[vi]
            ulod_flag = ulod_flags[vi]
            ulod_val = ulod_values[vi]
            vals = MaskedArray(dat * scale, mask=(dat == miss),
                               fill_value=miss)
            scale = scales[vi] = 1  # Set to 1 after applying
            tmpvar = self.variables[var] = PseudoNetCDFVariable(
                self, var, 'd', ('POINTS',), values=vals)
            tmpvar.units = unit
            tmpvar.standard_name = var
            tmpvar.missing_value = miss
            tmpvar.fill_value = miss
            tmpvar.scale = scale

            if hasattr(self, 'LLOD_FLAG'):
                tmpvar.llod_flag = llod_flag
                tmpvar.llod_value = llod_val

            if hasattr(self, 'ULOD_FLAG'):
                tmpvar.ulod_flag = ulod_flag
                tmpvar.ulod_value = ulod_val

        def dtime(s):
            return timedelta(seconds=int(s),
                             microseconds=(s - int(s)) * 1.E6)
        vtime = vectorize(dtime)
        tvar = self.variables[self.TFLAG]
        self._date_objs = (self._SDATE +
                           vtime(tvar).view(type=ndarray))


def ncf2ffi1001(f, outpath, mode='w', delim=', '):
    """
Arguments:
  f       - input file with 1-D variables and meta-data
  outpath - location to create output
  mode    - method for opening output file
  delim   - delimiter for data in output file
Returns:
  out     - output file (still open)
    """
    outfile = open(outpath, mode)
    check_for_attrs = ['PI_NAME', 'ORGANIZATION_NAME',
                       'SOURCE_DESCRIPTION', 'MISSION_NAME', 'VOLUME_INFO']
    missing_attrs = [k for k in check_for_attrs if k not in f.ncattrs()]
    if len(missing_attrs) > 0:
        warn('Missing import attributes filling with "Unknown": ' +
             ';'.join(missing_attrs))
    # header_keys = ("PI_CONTACT_INFO PLATFORM LOCATION ASSOCIATED_DATA " +
    #                "INSTRUMENT_INFO DATA_INFO UNCERTAINTY ULOD_FLAG " +
    #                "ULOD_VALUE LLOD_FLAG LLOD_VALUE DM_CONTACT_INFO " +
    #                "PROJECT_INFO STIPULATIONS_ON_USE OTHER_COMMENTS " +
    #                "REVISION").split()
    IGNORE_ATTRS = ['fmt', 'n_header_lines', 'PI_NAME', 'ORGANIZATION_NAME',
                    'SOURCE_DESCRIPTION', 'MISSION_NAME', 'VOLUME_INFO',
                    'SDATE', 'WDATE', 'TIME_INTERVAL', 'INDEPENDENT_VARIABLE',
                    'TFLAG']
    depvarkeys = [k for k in f.variables.keys() if k != f.INDEPENDENT_VARIABLE]
    myattrs = [k for k in f.ncattrs() if k not in IGNORE_ATTRS]
    print('%d, %d' % (len(myattrs) + len(depvarkeys) + 15, 1001), file=outfile)
    print(getattr(f, 'PI_NAME', 'Unknown'), file=outfile)
    print(getattr(f, 'ORGANIZATION_NAME', 'Unknown'), file=outfile)
    print(getattr(f, 'SOURCE_DESCRIPTION', 'Unknown'), file=outfile)
    print(getattr(f, 'MISSION_NAME', 'Unknown'), file=outfile)
    print(getattr(f, 'VOLUME_INFO', '1, 1'), file=outfile)
    print(f.SDATE, getattr(f, 'WDATE',
                           datetime.today().strftime('%Y, %m, %d')),
          file=outfile)
    print(getattr(f, 'TIME_INTERVAL', 0), file=outfile)
    print(f.INDEPENDENT_VARIABLE, file=outfile)
    print('%d' % len(depvarkeys), file=outfile)
    print(delim.join(['1' for k in depvarkeys]), file=outfile)
    print(delim.join([str(getattr(f.variables[k], 'missing_value', -999))
                      for k in depvarkeys]), file=outfile)
    for key, var in f.variables.items():
        if key == f.INDEPENDENT_VARIABLE:
            continue
        print(delim.join(
            [key, getattr(var, 'units', 'unknown')]), file=outfile)

    print(0, file=outfile)
    print(len(myattrs), file=outfile)
    for key in myattrs:
        print('%s: %s' % (key, getattr(f, key, '')), file=outfile)

    vals = [filled(f.variables[f.INDEPENDENT_VARIABLE][:]).ravel()]
    keys = [f.INDEPENDENT_VARIABLE]
    for key, var in f.variables.items():
        if key == f.INDEPENDENT_VARIABLE:
            continue
        keys.append(key)
        vals.append(filled(var[:]).ravel())

    print(delim.join(keys), file=outfile)
    for row in array(vals).T:
        row.tofile(outfile, format='%.6e', sep=delim)
        print('', file=outfile)

    return outfile


registerwriter('ffi1001', ncf2ffi1001)


class TestFfi1001(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF.testcase import icarttfiles_paths
        self.ffi1001path = icarttfiles_paths['ffi1001']

    def testNCF2FFI1001(self):
        import os
        ffi1001file = ffi1001(self.ffi1001path)
        outpath = self.ffi1001path + '.check'
        ncf2ffi1001(ffi1001file, outpath)
        newfile = ffi1001(outpath)
        for k, v in ffi1001file.variables.items():
            assert(k in newfile.variables)
            nv = newfile.variables[k]
            np.testing.assert_allclose(v[:], nv[:])

        os.remove(outpath)
