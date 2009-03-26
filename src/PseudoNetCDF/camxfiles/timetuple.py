__all__ = ['timerange', 'timediff', 'timeadd', 'cmp_time']

__doc__ = """
.. _timetuple
:mod:`timetuple` -- CAMx simple time functions
==============================================

.. module:: timetuple
   :platform: Unix, Windows
   :synopsis: Time tuple provided a simple interface for handling 
              CAMx time representations.  It has been deprecated
              and all reliance should be removed
.. moduleauthor:: Adam Hupp and Barron Henderson <barronh@unc.edu>
"""

from warnings import warn
warn(PendingDeprecationWarning("Time tuple will be replaced by datetime and dateutil"))

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

def timerange((date1,time1),(date2,time2), step=100, eod=2400.0):
    """Iterater of time tuples between start and end (not end inclusive)
    
    Should probably remove the != and use < if step + and > if step -
    """
    date1,time1=timeadd((date1,time1),(0,0),eod)
    date2,time2=timeadd((date2,time2),(0,0),eod)
    while (date1,time1)!=(date2,time2):
        yield date1,time1
        date1,time1 = timeadd((date1,time1),(0,step),eod)

def timediff((date1,time1),(date2,time2),  eod=2400.0):
    """Compares date tuples and returns difference in
    a time tuple
    """
    date3=date2-date1
    time3=time2-time1
    return date3*eod + time3

def timeadd((date1,time1),(date2,time2),eod=2400.0):
    """Adds time tuple 1 to time tuple 2.
    Time values greater than eod (dflt 2400)
    are converted to days
    
    Big question is hr 24 day 1
    """
    
    time1 = time1 + time2
    date1 = date1 + date2
    if time1 >= eod:
      time1 %= eod
      date1+=1
    if time1<0:
        time1 %= eod
        date1-=1
    return date1,time1

def cmp_time(lhs, rhs):
    """Compare two time tuples.
    Each arg is a 2-tuple, the first element being the hour of the day * 100,
    the second being the date
    """
    ltime, ldate = lhs
    rtime, rdate = rhs
    if ldate == rdate:
        return ltime - rtime
    else:
        return ldate - rdate
      
#
# Test cases
#        
import unittest

class CompareTime(unittest.TestCase):
    def testCompare(self):
        assert cmp_time((1.0, 2), (2.0,2)) < 0
        assert cmp_time((2.0, 2), (2.0,2)) == 0
        assert cmp_time((3.0, 2), (2.0,3)) < 0
        assert cmp_time((3.0, 3), (2.0,3)) > 0
