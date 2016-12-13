import numpy as np
def mda8(arr, axis = None, keepdims = True):
    """
    Daily-Maximum 8-hour average concentration
       - can be applied to any dimensions, but make sense with time
       - returns the day max of the a8
    arguments:
       arr
    """
    if hasattr(arr, 'mask'):
        arr_m = np.ma.masked_array(arr)
        nvalid = np.apply_along_axis(np.convolve, axis = axis, arr = ~np.ma.getarraymask(arr_m), v = [1]*8, mode = 'valid')
        arra8_m = np.ma.masked_where(nvalid < 6, np.apply_along_axis(np.convolve, axis = axis, arr = arr_m.filled(0), v = [1]*8, mode = 'valid')/nvalid)
        arrmda8 = np.ma.masked_invalid(daymax(arra8_m, axis = axis))
    else:
        arra8 = np.apply_along_axis(np.convolve, axis = axis, arr = arr, v = [1/8.]*8, mode = 'valid')
        arrmda8 = daymax(arra8, axis = axis)
    return arrmda8

def _dayfunc(func, arr, axis = None, keepdims = True):
    daystartend = np.arange(0, arr.shape[axis], 24)
    return func.reduceat(arr, daystartend, axis = axis)[:-1]

def daymax(arr, axis = None, keepdims = True):
    return _dayfunc(np.maximum, arr, axis = axis, keepdims = keepdims)

def daymin(arr, axis = None, keepdims = True):
    return _dayfunc(np.minimum, arr, axis = axis, keepdims = keepdims)

def daymean(arr, axis = None, keepdims = True):
    return _dayfunc(np.mean, arr, axis = axis, keepdims = keepdims)

def daystd(arr, axis = None, keepdims = True):
    return _dayfunc(np.std, arr, axis = axis, keepdims = keepdims)

def dayvar(arr, axis = None, keepdims = True):
    return _dayfunc(np.var, arr, axis = axis, keepdims = keepdims)

