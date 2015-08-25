from PseudoNetCDF.netcdf import NetCDFFile
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from PseudoNetCDF.MetaNetCDF import file_master
from warnings import warn

def mrgidx(ipr_paths, irr_paths, idx):
    if isinstance(irr_paths,str):
        irrf = NetCDFFile(irr_paths)
    else:
        irrf = file_master([NetCDFFile(irr_path) for irr_path in irr_paths])
    
    if isinstance(ipr_paths,str):
        iprf = NetCDFFile(ipr_paths)
    else:
        iprf = file_master([NetCDFFile(ipr_path) for ipr_path in ipr_paths])
        
    
    # Process and Reaction keys should exclude TFLAG
    pr_keys = [pr for pr in iprf.variables.keys() if pr not in ('TFLAG',)]
    rr_keys = [rr for rr in irrf.variables.keys() if rr not in ('TFLAG',)]
    
    # Attempt to order reactions by number
    # this is not necessary, but is nice and clean
    try:
        rr_keys = [(int(rr.split('_')[1]), rr) for rr in rr_keys]
        rr_keys.sort()
        rr_keys = [rr[1] for rr in rr_keys]
    except:
        warn("Cannot sort reaction keys")
    
    # Processes are predicated by a delimiter
    prcs = list(set(['_'.join(pr.split('_')[:-1]) for pr in pr_keys]))
    # Species are preceded by a delimiter
    spcs = list(set(['_'.join(pr.split('_')[-1:]) for pr in pr_keys]))
    
    # Select a dummy variable for extracting properties
    pr_tmp = iprf.variables[pr_keys[0]]
    
    # Create an empty file and decorate
    # it as necessary
    outf = PseudoNetCDFFile()
    outf.Species = "".join([spc.ljust(16) for spc in spcs])
    outf.Process = "".join([prc.ljust(16) for prc in prcs])
    outf.Reactions = "".join([rr_key.ljust(16) for rr_key in rr_keys])
    outf.createDimension("PROCESS", len(prcs))
    outf.createDimension("SPECIES", len(spcs))
    outf.createDimension("RXN", len(rr_keys))
    outf.createDimension("TSTEP", pr_tmp[:,0,0,0].shape[0])
    outf.createDimension("TSTEP_STAG", len(outf.dimensions["TSTEP"])+1)
    outf.createDimension("ROW", 1)
    outf.createDimension("LAY", 1)
    outf.createDimension("COL", 1)
    outf.createDimension("VAR", 3)
    outf.createDimension("DATE-TIME", 2)
    tflag = outf.createVariable("TFLAG", "i", ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag.__dict__.update(dict(units = "<YYYYJJJ,HHDDMM>", var_desc = 'TFLAG'.ljust(16), long_name = 'TFLAG'.ljust(16)))
    tflag[:,:,:] = iprf.variables['TFLAG'][:][:,[0],:]
    shape = outf.createVariable("SHAPE", "i", ("TSTEP", "LAY", "ROW", "COL"))
    shape.__dict__.update(dict(units = "ON/OFF", var_desc = "SHAPE".ljust(16), long_name = "SHAPE".ljust(16)))
    shape[:] = 1
    irr = outf.createVariable("IRR", "f", ("TSTEP", "RXN"))
    irr.__dict__.update(dict(units = pr_tmp.units, var_desc = "IRR".ljust(16), long_name = "IRR".ljust(16)))
    ipr = outf.createVariable("IPR", "f", ("TSTEP", "SPECIES", "PROCESS"))
    irr.__dict__.update(dict(units = pr_tmp.units, var_desc = "IPR".ljust(16), long_name = "IPR".ljust(16)))

    for rr, var in zip(rr_keys,irr.swapaxes(0,1)):
        var[:] = irrf.variables[rr][:][idx]
        
    for prc, prcvar in zip(prcs,ipr.swapaxes(0,2)):
        for spc, spcvar in zip(spcs,prcvar):
            try:
                spcvar[:] = iprf.variables['_'.join([prc,spc])][:][idx]
            except KeyError as es:
                warn(str(es))

    return outf
    
if __name__ == '__main__':
    from PseudoNetCDF.pncdump import pncdump
    from PseudoNetCDF.cmaqfiles.pa import mrgidx
    import sys
    
    x = ', '.join(sys.argv[1:])
    try:
        pncdump(eval('mrgidx(%s)' % (x,)), name = mrgidx)
    except:
        print >> sys.stderr, "\nUsage: python -m PseudoNetCDF.cmaqfiles.pa \"ipr_paths, irr_paths, idx\"\n  ipr_paths - a single string or list of strings\n  irr_paths - a single string or list of strings\n  idx - 4D numpy slice indexes a time series\n\n\nExample:\n  $ python -m PseudoNetCDF.cmaqfiles.pa \"['CMAQ_IPR1.nc','CMAQ_IPR2.nc'], 'CMAQ_IRR.nc', (slice(None),0,1,1)\""
