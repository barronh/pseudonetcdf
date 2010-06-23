__all__ = ['platform_is_bigendian', 'freadnumpy', 'freadstruct', 'fread', 'needs_byteswap', 'check_read', 'RecordFile', 'unpack_from_file', 'seek_to_record', 'read_into', 'writeline', 'OpenRecordFile', 'Int2Asc', 'Asc2Int']

__doc__ = """
.. _FortranFileUtil
:mod:`FortranFileUtil` -- Fortran File interaction functions
============================================================

.. module:: FortranFileUtil
   :platform: Unix, Windows
   :synopsis: Provides utilities for interacting with Fortran 
              "unformatted" binary files.  These files often 
              have records with integer buffers and formatted data.
              Often the data has a byte order different than the
              machine doing the reading. These utilities provide ways 
              for the user to get data and ignore the byte order and
              data buffers.
.. moduleauthor:: Adam Hupp
"""

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

##Distribution Packages
import unittest
import sys,struct,os.path
from operator import mul

##Site-Packages
#use numpy if available
#otherwise use numeric
from numpy import array,zeros,reshape,product,shape,transpose,fromfile

platform_is_bigendian = (sys.byteorder != 'little')

def freadnumpy(ifile,count,read_type,return_type=None,byteswap=False):
    if return_type is None:
        return_type = read_type

    if byteswap ^ platform_is_bigendian:
        endian = '>'
    else:
        endian = '<'
    rfmt= "%s%s" % (endian, return_type)
    fmt = "%s%s" % (endian, read_type)
    return array(fromfile(ifile,dtype=fmt,count=count),dtype=return_type)

def freadstruct(ifile, count, read_type, return_type=None, byteswap=False):
    if return_type is None:
        return_type = read_type

    if byteswap ^ platform_is_bigendian:
        endian = '>'
    else:
        endian = '<'

    fmt = "%s%d%s" % (endian, count, read_type)

    elem = unpack_from_file(fmt, ifile)
    
    return array(elem, return_type)

#numpy style fread is more than 2 x as fast as struct
#preserving freadstruct (previously fread) until further notice
fread=freadnumpy

def needs_byteswap(bigendian):
    """Determines if byteswapping is necessary for the given endianness
    For example, if we are reading bigendian data, call

    need_byteswap(True)

    to see if that data needs to be swapped on the current platform
    """

    return bigendian != platform_is_bigendian

def check_read(requested, result):
    """Checks that the result of a read from a file was of the expected size
    Raises EOFError if nothing was read, or IOError if the read was incomplete
    """
    if len(result) != requested and requested != 0:
        if result == '':
            raise EOFError()
        else:
            raise IOError("Incomplete read, requested %d recieved %d" %
                          (requested, len(result)))


class RecordFile(object):
    """
    This class provides an easy interface to treat unformated 
    Fortran files like like record files.
    
    TODO: count=-1 does not work in aread, off by one element?!?
    TODO: error checking in aread, prevent reading past end of record
    """
    def tell(self):
        """
        use file object tell
        """
        return self.infile.tell()
    
    def __init__(self, infile, bigendian=True):
        """
        Arguments:
        infile --  an open file-like object.  Must be random access.  Must be a real
          file if the RecordFile.aread method is going to be used.
        bigendian --  boolean, True if file is bigendian
        
        ***Assumes single 4 byte integer pad on either end
        """
        if type(infile)==str or type(infile)==unicode:
            self.infile = open(infile,'rb')
        elif type(infile)==file:
            self.infile = infile
        else:
            raise TypeError, "infile is of " + str(type(infile)) + ": and can only be string or file"
        
        # Is this platform a different byteorder than the file?
        self.byteswap =  needs_byteswap(bigendian)

        if bigendian:
            self.format_prefix = '>'
        else:
            self.format_prefix = '<'
        
        self.infile.seek(0, 2)
        self.length = self.infile.tell()

        self._newrecord(0)


    def seek(self, offset):
        if offset < self.length:
            posnow = self.infile.tell()
            relmov = offset - posnow
            self.infile.seek(relmov,1)
        else:
            raise ValueError, "File length = %i; record start requested %f" % (self.length,offset)
            
    def _newrecord(self, offset):
        """Move to a new record starting at the given offset
        """
        self.seek(offset)
        self.record_start = offset
        self.record_size = self.unpack("i")[0]
            

    def unpack(self, fmt):
        """Unpack a set of values, like the struct module
        """
        return unpack_from_file(self.format_prefix + fmt, self.infile)

    def read(self, fmt):
        """Unpack a set of values, like the struct module
        """
        result=unpack_from_file(self.format_prefix + fmt, self.infile)
        self.next()
        return result;

    def next(self):
        """Move to the next record.
        Returns True if sucessful, False if no more records available.
        """

        # 4 bytes for record size, start and end
        offset = self.record_start + self.record_size + 8 
        if offset < self.length:
            self._newrecord(offset)
            return True
        else:
            self.infile.seek(0, 2) # EOF
            return False
            
    def previous(self):
        """Move to the previous record.
        Returns True if sucessful, False if no more records available.
        """
        
        # 4 bytes for record size of this record
        # 4 bytes for record size of previous record
        if self.tell()==4:
            return False
        offset = (-8,-4)[self.eof()]
        
        # move to size byte of previous record
        self.infile.seek(offset,1)
        
        # Read size of previous record
        offset = self.unpack("i")[0]
        
        #actual offset = this  record start - size - previous record length + header
        offset = self.tell() - offset - 8
        if offset >= 0:
            self._newrecord(offset)
            return True
        else:
            self.infile.seek(0, 1) # EOF
            return False


    def aread(self, type, count):
        """Read a Numeric array

        Arguments:
        type -- a type code, one of those used by the struct or Numeric modules
        count -- number of elements to read.  If not specified,
                  defaults to all remaining bytes in current record
        """

        if count == -1:
            # NOTE: This does not work right!
            remain = self.record_size - (self.infile.tell() -
                                         self.record_start)
            count = remain/struct.calcsize(type)

        return fread(self.infile, count,
                     type, type,
                     self.byteswap)


    def eof(self):
        """Returns true if this RecordFile is at the end of the file
        """
        return self.infile.tell() == self.length

    def restart_record(self):
        """Move to beginning of record
        """
        self._newrecord(self.record_start)

def unpack_from_file(fmt, filein):
    """Like struct.unpack, but reads from a file instead of a string
    """
    fmtsize = struct.calcsize(fmt)
    data = filein.read(fmtsize)
    check_read(fmtsize, data)
    return struct.unpack(fmt, data)

def seek_to_record(rf,rid, fmt):
    """Searches for a record beginning with rid by unpacking
    the first struct.calcsize(fmt) bytes and comparing the
    results
    """
    while True:
        cid = rf.unpack(fmt)
        if rid == cid:
            rf.restart_record()
            break
        else:
            if not rf.next():
                raise ValueError, "Time %s not found" % str(rid)

def read_into(rf, dest, id_fmt, data_fmt='f'):
    """Read an array from a RecordFile into a Numeric array.
    I don't know how this will work if dest is anything but 2D.
    
    Arguments:
    rf -- a RecordFile instance
    dest -- a Numeric array, possibly a slice.

    There must be at least product(shape(dest)) elements
    left in the current record
    
    """

    if rf.eof(): raise EOFError()
    
    id = rf.unpack(id_fmt)
    rd = rf.aread(data_fmt,reduce(mul,dest.shape))
    dest[...] = rd.reshape(dest.shape)
    return id

def writeline(d,fmt,ForceBig=True):
    """writeline appends length integers
    and determines if byteswap is necessary
    """
    rlen=struct.calcsize(fmt)
    rfmt="i"+fmt+"i"
    if sys.byteorder=='little' or ForceBig:
        rfmt='>'+rfmt
    d=[i for i in d]
    d.insert(0,rlen)
    d.append(rlen)
    try:
        return struct.pack(rfmt,*d)
    except:
        print d
        raise

def OpenRecordFile(rf):
    """All CAMx files use FortranFileUtil.RecordFiles
    as inputs.  This function decides if the input was
    of the right type.  If not, it tries to make a RecordFile
    from the input.
    
    rf - str, unicode, file, RecordFile
    """
    if type(rf)==RecordFile:
        pass
    else:
        rf=RecordFile(rf)
    rf._newrecord(0)
    return rf

def Int2Asc(mspec):
    """Some CAMx input files have text stored as
    integers.  This function helps to undo that
    """
    spcname=""
    for c in mspec:
        spcname+= chr((((c-32)/256-32)/256-32)/256)
    return spcname

def Asc2Int(spcname):
    """Some CAMx output files need the text stored
    as integers.  This function helps to do that
    """
    mspec=[]
    for c in spcname:
        mspec.append(int((((((ord(c)*256)+32)*256)+32)*256)+32))
    return mspec

  
class TestFileUtils(unittest.TestCase):
    def setUp(self):
        from tempfile import TemporaryFile as tf
        self.tmpfile=tf(mode='w+b')
        
        #writing tempfile with Fortran unformatted binary format
        #1st line is 0-19 as floats
        #2nd line is 0-19 as ints
        #3rd line is 0-19 as strings
        self.tmpfile.write('\x00\x00\x00P'+'\x00\x00\x00\x00?\x80\x00\x00@\x00\x00\x00@@\x00\x00@\x80\x00\x00@\xa0\x00\x00@\xc0\x00\x00@\xe0\x00\x00A\x00\x00\x00A\x10\x00\x00A \x00\x00A0\x00\x00A@\x00\x00AP\x00\x00A`\x00\x00Ap\x00\x00A\x80\x00\x00A\x88\x00\x00A\x90\x00\x00A\x98\x00\x00'+'\x00\x00\x00P')
        self.tmpfile.write('\x00\x00\x00P'+'\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00\x04\x00\x00\x00\x05\x00\x00\x00\x06\x00\x00\x00\x07\x00\x00\x00\x08\x00\x00\x00\t\x00\x00\x00\n\x00\x00\x00\x0b\x00\x00\x00\x0c\x00\x00\x00\r\x00\x00\x00\x0e\x00\x00\x00\x0f\x00\x00\x00\x10\x00\x00\x00\x11\x00\x00\x00\x12\x00\x00\x00\x13'+'\x00\x00\x00P')
        self.tmpfile.write('\x00\x00\x00\x14'+'The quick brown fox '+'\x00\x00\x00\x14')
        self.tmprf=RecordFile(self.tmpfile)
        
    def testAdvance(self):
        self.tmprf._newrecord(0)
        self.assertEquals(self.tmprf.tell(),4)
        self.tmprf.next()
        self.assertEquals(self.tmprf.tell(),92)
        self.tmprf.next()
        self.assertEquals(self.tmprf.tell(),180)
        self.failIf(self.tmprf.next())
        self.assertEquals(self.tmprf.tell(),204)
        self.failIf(self.tmprf.next())
        self.tmprf.previous()
        self.assertEquals(self.tmprf.tell(),180)
        self.tmprf.previous()
        self.assertEquals(self.tmprf.tell(),92)
        self.tmprf.previous()
        self.assertEquals(self.tmprf.tell(),4)
        self.failIf(self.tmprf.previous())
        self.assertEquals(self.tmprf.tell(),4)
        
    def testFloat(self):
        from numpy import arange
        self.tmprf._newrecord(0)
        self.assert_((arange(20,dtype='f')==self.tmprf.aread('f',20)).all())

    def testReadInto(self):
        from numpy import arange,zeros
        dest=zeros((4,5),'f')
        
        self.tmprf._newrecord(0)

        #Necessary because written to anticipate fortran swapping of axes
        read_into(self.tmprf,dest,'','f')
        self.assert_((arange(20,dtype='f').reshape((4,5))==dest).all())
    
    def testInt(self):
        from numpy import arange
        self.tmprf._newrecord(0)
        self.tmprf.next()
        self.assert_((arange(20,dtype='i')==self.tmprf.aread('i',20)).all())
    
    def testSeek(self):
        self.tmprf._newrecord(0)
        seek_to_record(self.tmprf,(0.,1.,2.),"fff")
        self.assertEquals(self.tmprf.tell(),4)
        seek_to_record(self.tmprf,(0,1,2),"iii")
        self.assertEquals(self.tmprf.tell(),92)
        seek_to_record(self.tmprf,('The',),'3s')
        self.assertEquals(self.tmprf.tell(),180)
        seek_to_record(self.tmprf,('T','h','e',),'3c')
        self.assertEquals(self.tmprf.tell(),180)
        
        
    def testStr(self):
        from numpy import array
        self.tmprf._newrecord(0)
        self.tmprf.next()
        self.tmprf.next()
        self.assert_((array("The quick brown fox ")==self.tmprf.aread('S20',1)).any())
    
    def runTest(self):
        pass
if __name__ == '__main__':
  unittest.main()
