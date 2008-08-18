HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from numpy import ones,zeros,array
from Scientific.IO.NetCDF import NetCDFFile as ncf
from datetime import date,datetime

def toIoapi(obj,grid,outfile=None,vars=None,exclude=False,units='ppmV'):
    from Scientific.IO.NetCDF import NetCDFFile as ncf
    from datetime import datetime
    now=datetime.now
    
    if outfile==None:
        from tempfile import NamedTemporaryFile as ntf
        outfile=ntf()
    elif type(outfile)==str:
        outfile=file(outfile,'wb')
    elif type(outfile)!=file:
        raise ValueError, "outfile must be either string or file or none"
    if "w" not in outfile.mode:
        raise ValueError, "outfile must be writable"
    
    outfile=ncf(outfile.name,'w')
    #Configure variables to add
    if vars==None:
       vars=obj.spcnames
       
    if exclude:
       vars=[spc for spc in obj.spcnames if spc not in vars]

    #Add dimensions!
    dimensions=[
                 ('TSTEP', None),
                 ('DATE-TIME', 2),
                 ('LAY', obj.nlayers),
                 ('VAR', len(vars)),
                 ('ROW', obj.ny),
                 ('COL', obj.nx)
                ]

    #Configure timerange for output
    timerange=[(d,t) for d,t in obj.timerange()]
    newtimerange=[]
    for v in timerange:
        newtimerange.append([(v[0]+2000000,v[1]*10000/obj.time_step)]*len(vars))

    for k,v in dimensions:
        outfile.createDimension(k,v)
    
    #Add time variables
    t=outfile.createVariable('TFLAG','i',('TSTEP', 'VAR', 'DATE-TIME'))
    t.assignValue(newtimerange)
    setattr(t,'units','<YYYYDDD,HHMMSS>')
    setattr(t,'long_name',"TFLAG           ")
    setattr(t,'var_desc',"Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                ")
    #Add variables
    for v in vars:
       vo=outfile.createVariable(v.strip(),'f',('TSTEP','LAY','ROW','COL'))
       vo.assignValue(obj.variables[v].getValue())
       setattr(vo,'long_name',v.ljust(16))
       setattr(vo,'units',units.ljust(16))
       setattr(vo,'var_desc',v.ljust(80))
       outfile.sync()

    #Add global attributes
    spcnames=[spc.ljust(16) for spc in vars]
    now=now()
    cdate=now.year*1000+((now-datetime(now.year,1,1)).days+1)
    ctime=now.hour*10000+now.minute*100+now.second
    attributes=[
                ('IOAPI_VERSION', "$Id: @(#) ioapi library version 3.0 OpenMP enabled $                            " ),
                ('EXEC_ID', "????????????????                                                                " ),
                ('FTYPE', 1 ),
                ('CDATE', cdate ),
                ('CTIME', ctime ),
                ('WDATE', cdate ),
                ('WTIME', ctime ),
                ('SDATE', 2000000+obj.start_date ),
                ('STIME', int(obj.start_time*10000/obj.time_step) ),
                ('TSTEP', 10000 ),
                ('NTHIK', 1 ),
                ('NCOLS', obj.nx ),
                ('NROWS', obj.ny ),
                ('NLAYS', obj.nlayers ),
                ('NVARS', len(vars) ),
                ('GDTYP', 2 ),
                ('P_ALP', grid.alpha ),
                ('P_BET', grid.beta ),
                ('P_GAM', grid.gamma ),
                ('XCENT', grid.x_central ),
                ('YCENT', grid.y_central ),
                ('XORIG', grid.x_origin ),
                ('YORIG', grid.y_origin ),
                ('XCELL', grid.x_size ),
                ('YCELL', grid.y_size ),
                ('VGTYP', 5 ),
                ('VGTOP', numpy.array(grid.heights[-1],'f')),
                ('VGLVLS', grid.heights ),
                ('GDNAM', grid.name.ljust(16) ),
                ('UPNAM', "CAMxFiles       " ),
                ('VAR-LIST', ''.join(spcnames) ), # 16 chars
                ('FILEDESC', ("CAMxFile Avrg Conversion".ljust(80*60))),
                ('HISTORY', "" )
                ]
    for k,v in attributes:
        setattr(outfile,k,v)
        outfile.sync()
    return outfile

def MakeIoapi(nfile,dimensions,domain,layheights,gdnam,times,mode='w'):
    if type(nfile)==str:
        nfile=ncf(nfile,mode)
    for d,v in dimensions.iteritems():
        nfile.createDimension(d,v)
        nfile.sync()
    AddIoapiVariable(nfile,'TFLAG','<YYYYDDD,HHMMSS>',('TSTEP','VAR','DATE-TIME'),times,'i')
    nfile.sync()
    
    nfile.EXEC_ID="????????????????".ljust(80)
    nfile.IOAPI_VERSION = "$Id: @(#) ioapi library version 3.0 OpenMP enabled $".ljust(80)
    nfile.GDNAM = gdnam.ljust(16)
    nfile.UPNAM = "MassSelection".ljust(16)
    nfile.NCOLS=nfile.dimensions['COL']
    nfile.NROWS=nfile.dimensions['ROW']
    nfile.NLAYS=nfile.dimensions['LAY']
    nfile.NVARS=nfile.dimensions['VAR']
    nfile.FTYPE = 1 ;
    nfile.sync()

    nowd=date.today().year*1000+(date.today()-date(date.today().year,1,1)).days
    nowt=int(datetime.now().strftime("%H%M%S"))
    nfile.CDATE = getattr(nfile,'CDATE',nowd)
    nfile.CTIME = getattr(nfile,'CTIME',nowt)
    nfile.WDATE = nowd
    nfile.WTIME = nowt
    nfile.NTHIK = 1
    nfile.GDTYP = 2
    nfile.sync()

    nfile.P_ALP = domain.GetNormProjParm('standard_parallel_1')
    nfile.P_BET = domain.GetNormProjParm('standard_parallel_2')
    nfile.P_GAM = domain.GetNormProjParm("central_meridian")
    nfile.XCENT = domain.GetNormProjParm("central_meridian")
    nfile.YCENT = domain.GetNormProjParm("latitude_of_origin")
    nfile.XORIG = domain.GetNormProjParm("false_easting")
    nfile.YORIG = domain.GetNormProjParm("false_northing")
    nfile.XCELL = domain.GetLinearUnits()
    nfile.YCELL = domain.GetLinearUnits()
    nfile.sync()

    nfile.VGTYP = 6
    nfile.VGTOP = array(layheights ,'f')[-1]
    nfile.VGLVLS = array(layheights ,'f')
    nfile.FILEDESC = "Concentration file output".ljust(80)
    nfile.GDNAM = "".ljust(16)
    nfile.sync()

    setattr(nfile,'VAR-LIST',"")
    nfile.HISTORY = ""
    tflag=nfile.variables['TFLAG']
    nfile.TSTEP = tflag[1,0,1]-tflag[0,0,1]
    nfile.SDATE = tflag[0,0,0]
    nfile.STIME = tflag[0,0,1]
    nfile.sync()
    return nfile

def AddIoapiVariable(nfile,name,units,dimensions,value,typecode='f'):
    tempv=nfile.createVariable(name,typecode,dimensions)
    tempv.long_name=name
    tempv.var_desc=name
    tempv.units=units
    nfile.sync()
    tempv.assignValue(value)
    
    if name!='TFLAG':
        varlist=getattr(nfile,'VAR-LIST')
        varlist+=name.ljust(16)
        setattr(nfile,'VAR-LIST',varlist)
    nfile.sync()
    
