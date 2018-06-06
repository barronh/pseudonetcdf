__all__ = ['arl2time']


def arl2time(YY, MM, DD, HH, mm):
    """
    Convert ARL time pieces to a datetime object

    Parameters
    ----------
    YY : integer year (2 or 4 digits)
    MM : integer month (1-12)
    DD : integer day (1-31)
    HH : integer hour (0-23)
    mm : integer minute (0-59)

    Returns
    -------
    d : datetime object from datetime.strptime
    """
    from datetime import datetime
    strptime = datetime.strptime
    if YY < 100:
        yf = '%y'
    else:
        yf = '%Y'
    out = strptime('%02d-%02d-%02d %02d:%02d+0000' % (YY, MM, DD, HH, mm),
                   yf + '-%m-%d %H:%M%z')
    return out


def arl2timea(year, month, day, hour, minute):
    """
    Convert ARL time piece arrays to a datetime object array

    Parameters (all numpy array
    ----------
    YY : array integer year (2 or 4 digits)
    MM : integer month (1-12)
    DD : integer day (1-31)
    HH : integer hour (0-23)
    mm : integer minute (0-59)

    Returns
    -------
    d : array of datetime object from datetime.strptime
    """
    import numpy as np
    from datetime import datetime
    datei = (year.astype('l') * 100000000 +
             month.astype('l') * 1000000 +
             day.astype('l') * 10000 +
             hour.astype('l') * 100 +
             minute.astype('l'))
    datestrs = np.char.decode(np.char.add(datei.astype('S16'), b'+0000'))
    dates = np.array([datetime.strptime(d, '%y%m%d%H%M%z') for d in datestrs])
    return dates
