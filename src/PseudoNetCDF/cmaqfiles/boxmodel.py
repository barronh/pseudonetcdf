from ..sci_var import PseudoNetCDFFile
from numpy import array, dtype, fromstring, newaxis, vectorize
from datetime import timedelta, datetime

def box_model_conc(conc_path, start_datetime):
    lines = open(conc_path,'r').readlines()[8:]
    
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
        var.units = 'ppm'

    tflag = result.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    
    t_tflag = start_datetime+vectorize(lambda m: timedelta(minutes = m))(conc['Time'].astype('d'))
    tflag[:, :, 0] = vectorize(lambda d: int(d.strftime('%Y%j')))(t_tflag[:,newaxis])
    tflag[:, :, 1] = vectorize(lambda d: int(d.strftime('%H%M%S')))(t_tflag[:,newaxis])
    tflag.units = "<YYYYMMDD, HHMMSS>"
    tflag.long_name = 'TFLAG'.ljust(16)
    tflag.var_desc = 'TFLAG'.ljust(16)
 
    return result
 
def box_model_mrg(conc_path, irr_path, start_datetime):
    irr_values = open(irr_path,'r').readlines()
    irr_values = [[float(f) for f in l.split()] for l in irr_values]
    irr_values = array(irr_values)
    lines = open(conc_path,'r').readlines()[8:]
    
    spc_names = [spc_name.strip() for spc_name in lines[0][:-1].split('\t')]
    data_block = ''.join(lines[1:])
    
    
    line_format = dtype(dict(names=spc_names, formats=['f']*len(spc_names)))
    
    conc=fromstring(data_block, dtype='f', sep='\t').view(line_format)
    result = PseudoNetCDFFile()
    result.SDATE = int(start_datetime.strftime('%Y%j'))
    result.STIME = int(start_datetime.strftime('%H%M%S'))
    result.createDimension('TSTEP', conc.shape[0]-1)
    result.createDimension('PROCESS', 3)
    result.createDimension('REACTIONS', irr_values.shape[1]-2)
    result.createDimension('SPECIES', len(spc_names)-2)
    result.createDimension('VAR', 1)
    result.createDimension('DATE-TIME', 2)
    var = result.createVariable('IRR','f',('TSTEP', 'REACTIONS'))
    var.units = 'ppm'
    var[:,:] = irr_values[:,2:]
    
    var = result.createVariable('IPR','f',('TSTEP', 'SPECIES', 'PROCESS'))
    var.units = 'ppm'
    var[:,:,0] = conc.view('f').reshape(len(result.dimensions['TSTEP'])+1,len(result.dimensions['SPECIES'])+2)[:-1,1:-1]

    var[:,:,2] = conc.view('f').reshape(len(result.dimensions['TSTEP'])+1,len(result.dimensions['SPECIES'])+2)[1:,1:-1]
    var[:,:,1] = var[:,:,2] - var[:,:,0]

    tflag = result.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    
    t_tflag = start_datetime+vectorize(lambda m: timedelta(minutes = m))(conc['Time'].astype('d'))
    tflag[:, :, 0] = vectorize(lambda d: int(d.strftime('%Y%j')))(t_tflag[1:,newaxis])
    tflag[:, :, 1] = vectorize(lambda d: int(d.strftime('%H%M%S')))(t_tflag[1:,newaxis])
    result.Reactions = ''.join(['IRR_%d'.ljust(16) % ri for ri in range(1,len(result.dimensions['REACTIONS'])+1)])
    result.Species = ''.join([spc_name.ljust(16) for spc_name in spc_names if spc_name not in ['Time', 'AIR']])
    result.Process = 'INIT'.ljust(16)+'CHEM'.ljust(16)+'FCONC'.ljust(16)
    
    return result
