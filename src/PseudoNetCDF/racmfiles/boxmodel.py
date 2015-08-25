from ..sci_var import PseudoNetCDFFile
from numpy import array, dtype, fromstring, newaxis, vectorize
from datetime import timedelta, datetime

def box_model_irr(irr_path, start_datetime):
    retval = box_model_conc(irr_path, start_datetime)
    for spc, var in retval.variables.iteritems():
        if spc not in ('TFLAG',):
            var.units = '%s/h' % (var.units.strip(),)

    return retval
def box_model_conc(conc_path, start_datetime):
    lines = open(conc_path,'r').readlines()[3:]
    
    spc_names = [spc_name.strip() for spc_name in lines[0][:-1].split('\t')]
    data_block = ''.join(lines[1:])
    
    line_format = dtype(dict(names=spc_names, formats=['f']*len(spc_names)))
    conc=fromstring(data_block, dtype='f', sep='\t').view(line_format)
    
    result = PseudoNetCDFFile()
    result.SDATE = int(start_datetime.strftime('%Y%j'))
    result.STIME = int(start_datetime.strftime('%H%M%S'))
    result.createDimension('TSTEP', conc.shape[0])
    result.createDimension('VAR', 1)
    result.createDimension('DATE-TIME', 2)
    for spc in spc_names:
        var = result.createVariable(spc,'f',('TSTEP',))
        var[:] = conc[spc]
        var.long_name = spc.ljust(16)
        var.var_desc = spc.ljust(16)
        var.units = 'ppb'

    tflag = result.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    start_date = datetime.strptime(start_datetime.strftime('%Y%j'),'%Y%j')
    t_tflag = start_date+vectorize(lambda h: timedelta(hours = float(h)))(conc[spc_names[0]].astype('d'))
    tflag[:, :, 0] = vectorize(lambda d: int(d.strftime('%Y%j')))(t_tflag[:,newaxis])
    tflag[:, :, 1] = vectorize(lambda d: int(d.strftime('%H%M%S')))(t_tflag[:,newaxis])
    tflag.units = "<YYYYMMDD, HHMMSS>"
    tflag.long_name = 'TFLAG'.ljust(16)
    tflag.var_desc = 'TFLAG'.ljust(16)
 
    return result
 
def box_model_mrg(conc_path, irr_path, start_datetime):
    retval = PseudoNetCDFFile()
    conc = box_model_conc(conc_path, start_datetime)
    irr = box_model_irr(irr_path, start_datetime)
    sorted_rxns = [k for k in irr.variables.keys() if k not in ('TFLAG', 'TIME')]
    sorted_rxns.sort()
    sorted_spcs = [k for k in conc.variables.keys() if k not in ('TFLAG', 'TIME')]
    sorted_spcs.sort()

    nsteps = retval.createDimension('TSTEP', len(irr.dimensions['TSTEP'])-1)
    nrxns = retval.createDimension('REACTIONS', len(sorted_rxns))
    nrxns = len(retval.dimensions['REACTIONS'])
    nsteps = len(retval.dimensions['TSTEP'])
    retval.createDimension('SPECIES', len(sorted_spcs))
    retval.createDimension('PROCESS', 3)
    irr_vals = retval.createVariable('IRR','f',('TSTEP','REACTIONS'))
    irr_vals.units = 'ppb/h'
    irr_vals.long_name = 'IRR'.ljust(16)
    irr_vals.var_desc = 'IRR'.ljust(16)
    ipr_vals = retval.createVariable('IPR','f',('TSTEP','SPECIES', 'PROCESS'))
    ipr_vals.units = 'ppb/h'
    ipr_vals.long_name = 'IPR'.ljust(16)
    ipr_vals.var_desc = 'IPR'.ljust(16)

    for rxni, rxn_name in enumerate(sorted_rxns):
        vals = irr.variables[rxn_name][:]
        vals = vals.repeat(2,0)[1:-1].reshape(nsteps,2).mean(1)

        irr_vals[:,rxni] = vals

    irr_vals[0,:] = nan
    retval.Reactions = ''.join([k.ljust(16) for k in sorted_rxns])

    for spci, spc_name in enumerate(sorted_spcs):
        ipr_vals[:,spci,0]=conc.variables[spc_name][:-1]
        ipr_vals[:,spci,2]=conc.variables[spc_name][1:]
        ipr_vals[:,spci,1] = ipr_vals[:,spci,2]-ipr_vals[:,spci,0]

    retval.Species = ''.join([k.ljust(16) for k in sorted_spcs])
    retval.Process = 'INIT'.ljust(16)+'CHEM'.ljust(16)+'FCONC'.ljust(16)

    return retval
