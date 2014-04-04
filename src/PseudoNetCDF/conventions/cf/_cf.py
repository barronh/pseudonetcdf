import numpy as np
from datetime import datetime, timedelta
import re

def get_datetime_from_relative(var, format = '%Y-%m-%d %H:%M:%S'):
    units = var.units
    tunit, refstr = units.split(' since ')
    refstr = re.sub(r'\s*UTC', '', refstr)
    rdate = datetime.strptime(refstr, format)
    
    out = np.array([rdate + timedelta(**{tunit: float(dt)}) for dt in var[:]])
    return out
