import numpy as np
from PseudoNetCDF import PseudoNetCDFFile
class arlconcdump(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, path):
        try:
            cls(path, metaonly = True)
            return True
        except:
            return False
       
    def __init__(self, path, metaonly = False):
        self._infile = infile = open(path, 'rb')
        self._infile.seek(0,0)
        self._rec12345()
        if not metaonly: self._datarec()
    
    def _rec12345(self):
        """
        Record #1

        CHAR*4 Meteorological MODEL Identification
        INT*4 Meteorological file starting time (YEAR, MONTH, DAY, HOUR, FORECAST-HOUR)
        INT*4 NUMBER of starting locations
        INT*4 Concentration packing flag (0=no 1=yes) 
        """
        infile = self._infile
        rec1 = np.fromfile(infile, dtype = '>i,>S4,>i,>i,>i,>i,>i,>i,>i,>i', count = 1)[0]
        assert(rec1['f0'] == rec1['f9'])
        nloc = self.NSTARTLOCS = rec1['f7']
        pack = rec1['f8'] == 1
        self.METMODEL = rec1['f1'].decode()
        self.MET_YEAR = rec1['f2']
        self.MET_MONTH = rec1['f3']
        self.MET_DAY = rec1['f4']
        self.MET_HOUR = rec1['f5']
        self.MET_FORECAST_HOUR = rec1['f6']
        self.STARTING_LOCATIONS = rec1['f7']
        self.PACKED = rec1['f8']
        
        """
        Record #2 Loop to record: Number of starting locations

        INT*4 Release starting time (YEAR, MONTH, DAY, HOUR)
        REAL*4 Starting location and height (LATITUDE, LONGITUDE, METERS)
        INT*4 Release starting time (MINUTES) 
        """
        rec2 = np.fromfile(infile, dtype = '>i,>4i,>3f,>i,>i', count = nloc)
        assert((rec2['f0'] == rec2['f4']).all())
        self.createDimension('starts', nloc)
        sy = self.createVariable('START_YEAR', 'i', ('starts',), units = 'year')
        sy[:] = rec2['f1'][:, 0]
        sm = self.createVariable('START_MONTH', 'i', ('starts',), units = 'month')
        sm[:] = rec2['f1'][:, 1]
        sd = self.createVariable('START_DAY', 'i', ('starts',), units = 'day')
        sd[:] = rec2['f1'][:, 2]
        sh = self.createVariable('START_HOUR', 'i', ('starts',), units = 'hour')
        sh[:] = rec2['f1'][:, 3]
        slat = self.createVariable('START_LAT', 'f', ('starts',), units = 'degrees_north')
        slat[:] = rec2['f2'][:, 0]
        slon = self.createVariable('START_LON', 'f', ('starts',), units = 'degrees_east')
        slon[:] = rec2['f2'][:, 1]
        salt = self.createVariable('START_ALT', 'f', ('starts',), units = 'meters')
        salt[:] = rec2['f2'][:, 2]
        
        """
        Record #3

        INT*4 Number of (LATITUDE-POINTS, LONGITUDE-POINTS)
        REAL*4 Grid spacing (DELTA-LATITUDE,DELTA-LONGITUDE)
        REAL*4 Grid lower left corner (LATITUDE, LONGITUDE) 
        """
        rec3 = np.fromfile(infile, dtype = '>i,>2i,>2i,>2f,>i', count = 1)[0]
        assert(rec3['f0'] == rec3['f4'])
        
        nlats = self.NLATS = rec3['f1'][0]
        nlons = self.NLONS = rec3['f1'][1]
        self.DELTA_LAT = rec3['f2'][0]
        self.DELTA_LON = rec3['f2'][1]
        self.LLCRNR_LAT = rec3['f3'][0]
        self.LLCRNR_LON = rec3['f3'][1]
        """
        Record #4

        INT*4 NUMBER of vertical levels in concentration grid
        INT*4 HEIGHT of each level (meters above ground) 
        """
        tmp = np.fromfile(infile, dtype = '>i,>i', count = 1)[0]
        nlays = self.NLAYS = tmp['f1']
        infile.seek(-8, 1)
        rec4 = np.fromfile(infile, dtype = '>i,>i,>{}i,>i'.format(nlays), count = 1)[0]
        assert(rec4['f0'] == rec4['f3'])
        self.createDimension('layer', nlays)
        var = self.createVariable('layer', 'i', ('layer',))
        var.units = 'meters agl'
        var[:] = rec4['f2']
        
        """
        Record #5

        INT*4 NUMBER of different pollutants in grid
        CHAR*4 Identification STRING for each pollutant 
        """
        tmp = np.fromfile(infile, dtype = '>i,>i', count = 1)
        npols = self.NPOLS = tmp['f1'][0]
        infile.seek(-8, 1)
        rec5 = np.fromfile(infile, dtype = '>i,>i,>({},)S4,>i'.format(npols), count = 1).squeeze()
        assert(rec5['f0'] == rec5['f3'])
        self.POLLUTANTS = b','.join(rec5['f2']).decode()

    def _datarec(self):
        infile = self._infile
        nlats = self.NLATS
        nlons = self.NLONS
        nlays = self.NLAYS
        npols = self.NPOLS
        pack = self.PACKED == 1
        now = infile.tell()
        infile.seek(0, 2)
        total = infile.tell()
        infile.seek(now, 0)
        datasize = total - now
        datasize
        thdr = np.dtype('>i,>6i,>i')
        nopackfmt = np.dtype('>i,>S4,>i,>({},{})f,>i'.format(nlats, nlons))
        outs = []
        pols = []
        lays = []
        starts = []
        stops = []
        ctmp = np.zeros((nlats, nlons), dtype = 'f')
        while infile.tell() != total:
            """
            Record #6 Loop to record: Number of output times

            INT*4 Sample start (YEAR MONTH DAY HOUR MINUTE FORECAST) 
            """
            start = np.fromfile(infile, dtype = thdr, count = 1)
            starts.append(start)
            """
            Record #7 Loop to record: Number of output times

            INT*4 Sample stop (YEAR MONTH DAY HOUR MINUTE FORECAST) 
            """
            stop = np.fromfile(infile, dtype = thdr, count = 1)
            stops.append(stop)
            """
            Record #8 Loop to record: Number levels, Number of pollutant types

            CHAR*4 Pollutant type identification STRING
            INT*4 Output LEVEL (meters) of this record

            No Packing (all elements)
            REAL*4 Concentration output ARRAY

            Packing (only non-zero elements)

            INT*4 Loop non-zero elements
            INT*2 First (I) index value
            INT*2 - Second (J) index value
            REAL*4 - Concentration at (I,J)**_ 
            """
            for i in range(npols * nlays):
                if pack:
                    forfmt = np.dtype('>i,>S4,>i,>i')
                    tmp = np.fromfile(infile, forfmt, count = 1)
                    myb = tmp['f0'][0]
                    myp = tmp['f1'][0]
                    myl = tmp['f2'][0]
                    pols.append(myp)
                    lays.append(myl)
                    myc = tmp['f3'][0]
                    infile.seek(-16, 1)
                    c = np.zeros((nlats, nlons), dtype = 'f')
                    tmp = np.fromfile(infile, np.dtype([('f0', forfmt, 1), ('data', np.dtype(dict(names = 'IJC', formats = '>i2,>i2,>f'.split(','))), myc), ('f1', '>i',1)]), count = 1)
                    J = tmp[0]['data']['J']
                    I = tmp[0]['data']['I']
                    C = tmp[0]['data']['C']
                    c[J, I] = C
                else:
                    c = np.fromfile(infile, dtype = nopackfmt, count = 1)[0]
                outs.append(c)
        datablock = np.array(outs).reshape(-1, nlays, npols, nlats, nlons)
        ntimes = datablock.shape[0]
        pols = np.array(pols).reshape(ntimes, nlays, npols)
        lays = np.array(lays).reshape(ntimes, nlays, npols)
        self.createDimension('time', ntimes)
        self.createDimension('layer', nlays)
        self.createDimension('latitude', nlats)
        self.createDimension('longitude', nlons)
        assert((lays[:, :, [0]] == lays[:, :, :]).all())
        assert((pols[:, [0], :] == pols[:, :, :]).all())
        for pi, pol in enumerate(pols[0,0]):
            var = self.createVariable(pol.decode(), 'f', ('time', 'layer', 'latitude', 'longitude'))
            var.units = 'arbitrary'
            var.description = pol.decode()
            var[:] = datablock[:, :, pi]
if __name__ == '__main__':
    print(arlconcdump.isMine('cdump24'))
    f = arlconcdump('cdump24')
