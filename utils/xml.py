HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

def validate_xml(xml_filename, dtd_filename=None,Raiser=True):
    """Validate a given XML file with a given external DTD.

    If the XML file is not valid, an error message will be printed
    to sys.stderr, and the program will be terminated with a non-zero
    exit code.  If the XML file is valid, nothing will be printed.
    """
    from xml.parsers.xmlproc import utils, xmlval, xmldtd,xmlproc
    

    if dtd_filename!=None:
        parser=xmlproc.XMLProcessor()
        try:
            xmldtd.load_dtd(dtd_filename)
            parser.dtd = parser.dtd = parser.ent = dtd
            parser.set_application(xmlval.ValidatingApp(dtd, parser))
        except:
            parser=xmlval.XMLValidator()
            raise UserWarning, "DTD was unavailable"
    else:
        parser=xmlval.XMLValidator()
        
    if Raiser:
        parser.set_error_handler(utils.ErrorRaiser(parser))
    else:
        parser.set_error_handler(utils.ErrorPrinter(parser))
    try:
        parser.parse_resource(xml_filename)   
    except:
        raise UserWarning, "XML was not validated"

def SingleNode(path,node):
    from xml.xpath import Evaluate
    return Evaluate(path,node)[0]

def SingleNodeTextValue(path, node):
    from xml.xpath import Evaluate
    if path[-7:] != "/text()":
        path += "/text()"
    return SingleNode(path,node).nodeValue
    
def job_info_xml(jobfile,dtidx=0):
    """Reader function for the pyPA xml format
    
    This returns an object with properties of 
    the xml file for easy access
    """
    from xml.xpath import Evaluate
    from xml.dom import minidom
    from pyPA.utils.util import validate_xml,AttrDict
    import sys
    #Check the xml against it's dtd
    try:
        validate_xml(jobfile)
    except:
        print >> sys.stderr, "XML Not validated"
        
    #open the xml file and get the documentElement
    kvdoc = minidom.parse(open(jobfile)).documentElement
    
    #create an object to store attributes
    result = AttrDict()
    
    #find the kvfile text node in the document element
    result.kvfile = SingleNodeTextValue("kvfile", kvdoc)

    #find the zpfile text node in the document element
    result.zpfile = SingleNodeTextValue("zpfile", kvdoc)

    #find the meteorology case
    result.meteorology = SingleNodeTextValue("meteorology", kvdoc)
    
    #find the emissions case
    result.emissions = SingleNodeTextValue("emissions", kvdoc)

    #find the zpfile text node in the document element
    try:
        result.extfile = SingleNodeTextValue("extfile", kvdoc)
    except:
        result.extfile = 'ext_pa.'+result.meteorology+'.'+result.emissions+'.'+SingleNodeTextValue('focus_domain/name',kvdoc)+'.%s.ext'
    
    #find the zpfile text node in the document element
    try:
        result.mrgfile = SingleNodeTextValue("mrgfile", kvdoc)
    except:
        result.mrgfile = 'mrgaloft.'+result.meteorology+'.'+result.emissions+'.'+SingleNodeTextValue('focus_domain/name',kvdoc)+'.%s.ext'
    
    
    #find the model name and version
    model=SingleNode("model",kvdoc)
    result.model=AttrDict()
    result.model.name=SingleNodeTextValue("name",model)
    result.model.version=SingleNodeTextValue("version",model)
    
    #find the kvfile text node in the document element
    result.ptdelim = SingleNodeTextValue("pointdelim", kvdoc)

    #find the kvfile text node in the document element
    result.coordelim = SingleNodeTextValue("coordelim", kvdoc)
    
    #find the irrfile text node in the document element
    try:
      result.irr_file = SingleNodeTextValue("irrfile", kvdoc)
    except:
      print "No irr file provided"
    
    #find the irrfile text node in the document element
    try:
      result.ipr_file = SingleNodeTextValue("iprfile", kvdoc)
    except:
      print "No ipr file provided"
    
    
    date=Evaluate("focus_domain/date",kvdoc)[dtidx]
    #start and end hr and dt should be specified in the file
    result.stime=(
                    int(SingleNodeTextValue("startdt", date)),
                    float(SingleNodeTextValue("starthr", date))
                   )
    result.etime=(
                    int(SingleNodeTextValue("enddt", date)),
                    float(SingleNodeTextValue("endhr", date))
                   )
    result.datestr=SingleNodeTextValue("datestr",date)

    #the focus domain is a polygon or line 
    #described by a list of points.  The coordinates
    #are 1-based like the domain, but the array is 0 based
    #so we subtract one from the coordinates
    try:
        result.ijpbl=dict([
                                 (
                                    SingleNodeTextValue("id",i),
                                      [add_tuple(toints(j,result.coordelim),(-1,-1,0)) 
                                          for j in SingleNodeTextValue("ijpbl", i).split(result.ptdelim)]
                                  ) 
                                      for i in Evaluate("hour",date)]
                            )

    except IndexError:
        result.geometry=dict(
                              [
                                  (
                                      SingleNodeTextValue("id", i),
                                      [add_tuple(toints(j,result.coordelim),(-1,-1)) 
                                          for j in SingleNodeTextValue("geometry", i).split(result.ptdelim)]
                                  ) 
                                      for i in Evaluate("hour",date)]
                          )
        try:
            result.pbl=dict(
                                  [
                                      (
                                          SingleNodeTextValue("id", i),
                                          SingleNodeTextValue("pbl", i)
                                      ) 
                                          for i in Evaluate("hour",date)]
                              )        
        except:
            print >> sys.stderr, "PBL not provided"
    pagrid=SingleNode("pagrid",kvdoc)
    pagrido=grid(
        int(SingleNodeTextValue("xorg",pagrid)),
        int(SingleNodeTextValue("ncols",pagrid)),
        int(SingleNodeTextValue("yorg",pagrid)),
        int(SingleNodeTextValue("nrows",pagrid)),
        int(SingleNodeTextValue("size",pagrid)),
        int(SingleNodeTextValue("layers",pagrid))
        )
    result.volume_shape=(pagrido.ncols,pagrido.nrows,pagrido.nlays)
    try:
        kvgrid=SingleNode("kvgrid",kvdoc)
        kvgrido=grid(
            int(SingleNodeTextValue("xorg",kvgrid)),
            int(SingleNodeTextValue("ncols",kvgrid)),
            int(SingleNodeTextValue("yorg",kvgrid)),
            int(SingleNodeTextValue("nrows",kvgrid)),
            int(SingleNodeTextValue("size",kvgrid)),
            int(SingleNodeTextValue("layers",kvgrid))
            )
        
        #Geometry is 0 indexed, but grid objects are 1 indexed
        #this lambda 1 indexes the request; converts and then 0 indexes the return
        result.kvgridij=lambda i,j: add_tuple(gridconverter(pagrido,kvgrido).ij(*add_tuple((i,j),(1,1))),(-1,-1))
        result.kv_shape=(kvgrido.ncols,kvgrido.nrows,kvgrido.nlays)

    except:
        result.kvgridij=lambda i,j: (i,j)
        result.kv_shape=result.volume_shape

    for kf in [k for k in result.keys() if k[-4:]=='file']:
        result[kf]=result[kf] % result.datestr
    #The krange is now being set automatically
    #do we need the ability to set this value
    result.krange = (0,result.volume_shape[2])

    return result
