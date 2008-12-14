HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

__all__ = ['irr2ncf', 'ipr2ncf']

def irr2ncf(irrinf,irroutf,lay=None,row=None,col=None,verbose=False):
    """
    Convert the CAMx irr files for the 4-k and 1-k simulations to the NetCDF format
    """
    if type(irrinf)==str:
        try:
            irrinf=irr_memmap(irrinf)
        except OverflowError:
            irrinf=irr(irrinf)

    if type(irroutf)==str:
        irroutf=ncf(irroutf,'w+')

    if lay==None:
        lay=int(irrinf.dimensions['LAY'])
    
    if row==None:
        row=int(irrinf.dimensions['ROW'])
    
    if col==None:
        col=int(irrinf.dimensions['COL'])
    
    if verbose:
        print >> sys.stderr, "Add irr Globals"
    irroutf.UPNAM = "irr FILE         " ;
    irroutf.GDNAM = "" ;
    irroutf.IOAPI_VERSION = "2.2 2003141 (May 21, 2003)" ;
    irroutf.FILEDESC = irrinf.rffile.infile.name
    irroutf.VGTYP = 2 ;
    irroutf.GDTYP = 2 ;
    irroutf.FTYPE = 1 ;
    irroutf.HISTORY = "New File written by irr2ncf" ;
    irroutf.EXEC_ID = "" ;
    irroutf.NLAYS = lay ;
    irroutf.NROWS = row ;
    irroutf.NCOLS = col ;
    irroutf.NVARS = 2 ;
    irroutf.NTHIK = 1 ;
    try:
        irroutf.SDATE = irrinf.SDATE ;
        irroutf.STIME = irrinf.STIME ;
        irroutf.TSTEP = irrinf.TSTEP ;
    except:
        irroutf.SDATE = irrinf.start_date ;
        irroutf.STIME = irrinf.start_time ;
        irroutf.TSTEP = irrinf.time_step ;
        

    if verbose:
        print >> sys.stderr, "Add irr Dimensions"
        
    irroutf.createDimension('VAR',int(irroutf.NVARS))
    irroutf.createDimension('LAY',int(irroutf.NLAYS))
    irroutf.createDimension('ROW',int(irroutf.NROWS))
    irroutf.createDimension('COL',int(irroutf.NCOLS))
    irroutf.createDimension('TSTEP',None)

    if verbose:
        print >> sys.stderr, "Create Variables"
        
    rxn_iter=[i for i in irrinf.variables.keys()]
    
    try:
        j=array(irrinf.memmaps['J'])
        i=array(irrinf.memmaps['I'])
        rowwindow=slice(j.min()-1,j.max())
        colwindow=slice(i.min()-1,i.max())
    except:
        dom=irrinf.activedomain
        rowwindow=slice(dom['jstart']-1,dom['jend'])
        colwindow=slice(dom['istart']-1,dom['iend'])

        
    for rxn in rxn_iter:
            print >> sys.stderr, rxn
            tmp_var=irroutf.createVariable(rxn,'f',('TSTEP','LAY','ROW','COL'))
            tmp_var.long_name=rxn
            tmp_var.var_desc=tmp_var.long_name
            tmp_var.units='ppm'
            irroutf.sync()

    if verbose:
        print >> sys.stderr, "irr Variables Fill"

    for rxn in rxn_iter:
            print >> sys.stderr, rxn
            tmp_var=irroutf.variables[rxn]
            tmp_var.assignValue(0)
            tmp_var[:,:,rowwindow,colwindow]=irrinf.variables[rxn]
            irroutf.sync()
    
    irrinf.close()
    if verbose:
        print >> sys.stderr, "irr Done"
        
    return irroutf
    
def ipr2ncf(iprinf,iproutf,lay=None,row=None,col=None,verbose=False):
    """
    Convert the CAMx ipr files for the 4-k and 1-k simulations to the NetCDF format
    """
    if type(iprinf)==str:
        try:
            iprinf=ipr_memmap(iprinf)
        except OverflowError:
            iprinf=ipr(iprinf)

    if type(iproutf)==str:
        iproutf=ncf(iproutf,'w+')

    if lay==None:
        lay=int(iprinf.dimensions['LAY'])
    
    if row==None:
        row=int(iprinf.dimensions['ROW'])
    
    if col==None:
        col=int(iprinf.dimensions['COL'])
    
    if verbose:
        print >> sys.stderr, "Add IPR Globals"
    iproutf.UPNAM = "IPR FILE         " ;
    iproutf.GDNAM = "" ;
    iproutf.IOAPI_VERSION = "2.2 2003141 (May 21, 2003)" ;
    iproutf.FILEDESC = iprinf.rffile.infile.name
    iproutf.VGTYP = 2 ;
    iproutf.GDTYP = 2 ;
    iproutf.FTYPE = 1 ;
    iproutf.HISTORY = "New File written by ipr2ncf" ;
    iproutf.EXEC_ID = "" ;
    iproutf.NLAYS = lay ;
    iproutf.NROWS = row ;
    iproutf.NCOLS = col ;
    iproutf.NVARS = 2 ;
    iproutf.NTHIK = 1 ;
    try:
        iproutf.SDATE = iprinf.SDATE ;
        iproutf.STIME = iprinf.STIME ;
        iproutf.TSTEP = iprinf.TSTEP ;
    except:
        iproutf.SDATE = iprinf.start_date ;
        iproutf.STIME = iprinf.start_time ;
        iproutf.TSTEP = iprinf.time_step ;

    if verbose:
        print >> sys.stderr, "Add IPR Dimensions"
        
    iproutf.createDimension('VAR',int(iproutf.NVARS))
    iproutf.createDimension('LAY',int(iproutf.NLAYS))
    iproutf.createDimension('ROW',int(iproutf.NROWS))
    iproutf.createDimension('COL',int(iproutf.NCOLS))
    iproutf.createDimension('TSTEP',None)

    if verbose:
        print >> sys.stderr, "Create Variables"
        
    try:
        prc_iter=[i for i in iprinf.proc_dict.keys()]
    except:
        prc_iter=[i for i in iprinf.ipr_record_type.names]
    
    try:
        pi=prc_iter.index('SPC')
        prc_iter.pop(pi)
    except ValueError:
        pass
    
    try:
        j=array(irrinf.memmaps['J'])
        i=array(irrinf.memmaps['I'])
        rowwindow=slice(j.min()-1,j.max())
        colwindow=slice(i.min()-1,i.max())
    except:
        dom=iprinf.activedomain
        rowwindow=slice(dom['jstart']-1,dom['jend'])
        colwindow=slice(dom['istart']-1,dom['iend'])

        
    for spc in iprinf.spcnames:
        for prc in prc_iter:
            prc_spc=prc.strip()+'_'+spc.strip()
            print >> sys.stderr, prc_spc
            tmp_var=iproutf.createVariable(prc_spc,'f',('TSTEP','LAY','ROW','COL'))
            tmp_var.long_name=prc.strip()+" and "+spc.strip()
            tmp_var.var_desc=tmp_var.long_name
            tmp_var.units='umol/m3'
            iproutf.sync()

    if verbose:
        print >> sys.stderr, "IPR Variables Fill"

    for spc in iprinf.spcnames:
        for prc in prc_iter:
            prc_spc=prc.strip()+'_'+spc.strip()
            print >> sys.stderr, prc_spc
            tmp_var=iproutf.variables[prc_spc]
            tmp_var.assignValue(0)
            tmp_var[:,:,rowwindow,colwindow]=iprinf.variables[prc_spc]
            iproutf.sync()
    
    iprinf.close()
    if verbose:
        print >> sys.stderr, "IPR Done"
        
    return iproutf

