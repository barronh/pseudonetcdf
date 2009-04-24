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
    lines = file(conc_path,'r').readlines()[3:]
    
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
 
