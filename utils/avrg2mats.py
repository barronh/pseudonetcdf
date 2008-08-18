HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

"""Create model prediction "MATS-format" ascii file

usage:
from avrg2mats import *

This file official source is $HeadURL$
"""

ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
Rnum = RevisionNum.split()[1]
ChangedBy  = "$LastChangedBy$"

# import libraries
# import standard modules
import csv

# site-packages
import gdal
from mx.DateTime import DateTimeDeltaFrom, DateTimeFrom, oneDay, oneHour, ISO
import numpy

# pyPASS utility module 
import DomainGrid
from XML2Dicts import GetXMLRootAttrs, XML2DomainDict

# import pyPA utils
#from pyPA.utils.CAMxFiles import gridded_emissions as ge
from CAMxFiles import gridded_emissions as ge

all = [ "getDataArray", "writeMATSmodel" ] 

def getDataArray(
        inputtemplate="test/camx420_cb4.%d.hgb8h.b1b.psito2n2.TCEQuh1_eta_tke.PA.favrg.hgbpa_04k", 
        #daterange=range(20000818,20000832)+range(20000901,20000907)
        daterange=[20000818]
                 ):
    
    e = ge(inputtemplate % daterange[0])
    HH          = 24 * len(daterange)
    arrayshape  = ( HH      ,
                    e.nx    , # cols
                    e.ny      # rows
                  )
    fullo3array = numpy.empty(arrayshape, dtype="float32")
    o38hdata    = numpy.empty(arrayshape, dtype="float32")
    
    print "Extracting 1HO3 concs from avrg files..."
    for i, date in enumerate(daterange) :
        
        print "  file:", inputtemplate % date
        f = ge(inputtemplate % date)
        
        o3indx = f.spcnames.index("O3".ljust(10))
        o3data = f.getArray(nspec=slice(o3indx,o3indx+1),krange=slice(1,2))
        
        fullo3array[i*24:(i+1)*24] = o3data[0:24,0,:,:,0]
    
    print "Calculating 8HO3 average concs..."
    for h in range(HH-6) :
        o38hdata[h] = numpy.average(fullo3array[h:h+8],0)
    
    max8ho3data = numpy.empty((len(daterange),e.nx,e.ny),"float32")
    print "Calculating maximum 8HO3 concs..."
    for i, date in enumerate(daterange) :
        print " ", date
        max8ho3data[i] = numpy.amax(o38hdata[i*24:(i+1)*24],0)
        #print numpy.amax(o38hdata[i*24:(i+1)*24],0).max(), max8ho3data[i].max()
    
    print "8HO3 data array shape:", max8ho3data.shape
    return max8ho3data
    
def getLatLongArray(
                   domainfile="dom.xml"                 ,
                   domain=0                             ,
                   ):
    
    # get domain information
    dom = [ d for d in XML2DomainDict(domainfile) ][domain]
    print "Retrieved information for domain %d:", domain
    print "   %s" % dom['identity.name']
    print "   Horizontal grid resolution:", dom['hgrid.resolution']
    
    # get domain file attributes
    domattrs = GetXMLRootAttrs(domainfile)
    
    # get domain grid information (including origins)
    dominfo     = DomainGrid.DomainGrid(
                                         dom                    ,
                                         domattrs['gc_prj4_str'],
                                         domattrs['pj_prj4_str']
                                        )
    # get origin info
    r  = dominfo.resolution            * 1000  # km to m
    eo = (dominfo.orig_easting  + 0.5) * 1000  # km to m
    no = (dominfo.orig_northing + 0.5) * 1000  # km to m
    
    # get geographic coordinate object
    gc = gdal.osr.SpatialReference()
    gc.ImportFromProj4(domattrs['gc_prj4_str'])
    
    # get projection coordinate object
    pj = gdal.osr.SpatialReference()
    pj.ImportFromProj4(domattrs['pj_prj4_str'])
    pj.CopyGeogCSFrom(gc)
    
    # create transformer object LCC --> lat long
    LCCtoLatLon = gdal.osr.CoordinateTransformation(pj, gc)
    transformer = LCCtoLatLon.TransformPoints
    
    # create list of lcc points
    lccpoints  = [(eo+r*i,no+r*j) for i in range(dominfo.num_columns) for j in range(dominfo.num_rows)]
    latlonvals = transformer(lccpoints)
    
    # create array with shape nx, ny, 3=lat/lon/ele
    latlongarray = numpy.reshape(latlonvals, (dominfo.num_columns,dominfo.num_rows,3))
    
    return latlongarray
    
def writeMATSmodel(
                   max8ho3data=numpy.zeros((1,1,1),"float32") , 
                   latlongarray=numpy.zeros((1,1,3), "float32"),
                   filename="test.csv"                  , 
                   firstdate="2000-08-18"               ,
                   ):
    
    (ndays, nx, ny) = max8ho3data.shape
    firstday = DateTimeFrom(firstdate)
    
    print "opening csv output file:", filename
    f = open(filename, "w")
    w = csv.writer(f, doublequote=False, quoting=csv.QUOTE_NONE, escapechar='"')
    
    print "writing header information"
    w.writerow(["Day"])
    #headerline = "  _ID,_TYPE,   LAT ,  LONG  ,    DATE  , O3"
    headsize   = [ 5, 5, 7, 8, 10, 3 ]
    headers    = ["_ID","_TYPE","LAT ","LONG  ","DATE  ","O3"]
    w.writerow( [ entry.rjust(headsize[i]) for i, entry in enumerate(headers) ] )
    f.flush()
    
    #outputarray = resize(array("", dtype="|S11"),max8ho3data.shape) # 7 len str
    #final output row should look like this:
    #  4004,"",28.800378, -99.133951,20010501, 49.0873
    datasize = [ 6, 1, 8, 11, 8, 8 ]
    
    # lat & long data: center of cell
    #lat = "fixme" ; long= "metoo"
    
    max8ho3data = numpy.where(max8ho3data <= 0.0, -9, max8ho3data*1000)
    
    print "Writing data to file"
    for col in range(nx) :
        for row in range(ny) :
            for day in range(len(max8ho3data)) :
                content = [
     (col+1)*1000 + row+1                       ,
    '"'                                         ,
    latlongarray[col,row,1]                     ,
    latlongarray[col,row,0]                     ,
    (firstday + oneDay * day).strftime("%Y%m%d"),
    "%6.4f" % max8ho3data[day,col,row]          ,
                          ]
                for i, val in enumerate(content) :
                    content[i] = str(val)[0:datasize[i]+1].rjust(datasize[i])
                if int(float(content[-1].strip())) == -9 : content[-1] = "-9"
                w.writerow( content )
                f.flush()
    f.close()