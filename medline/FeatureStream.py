"""A class for rapid iteration over the records in Medline."""

from bsddb import db
import logging
import numpy as nx
import struct


                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""

        

def Date2Integer(date):
    """Given (year,month,day), return the integer representation"""
    return date[0]*10000 + date[1]*100 + date[2]


def Integer2Date(intdate):
    """Given a YYYYMMDD integer, return (year,month,day)"""
    return intdate//10000, (intdate%10000)//100, intdate%100


class PlainFeatureStream:
    """Class for reading/writing a binary stream of Medline records, consisting
    of PubMed ID, record completion date and a vector of Feature IDs for
    features present in the record.  This stream is 
    
    @ivar stream: File object underlying the stream.
    
    @ivar filename: Path to the file holding the feature stream.
    
    @ivar ftype: Numpy integer type for representing features.
    
    @ivar rdonly: Boolean for whether to open the stream read-only or appendable.
    """
    
    
    max_bytes = 2000
    """Max feature vector bytes to read (more implies a stream error)."""


    def __init__(self, filename, ftype, rdonly):
        """Initialise the stream."""
        self.filename = filename
        self.ftype = ftype
        self.rdonly = rdonly
        if not filename.exists():
            filename.touch()
        self.stream = open(filename, "rb" if rdonly else "rb+")
        self.stream.seek(0,2) # Go to end of file


    def close(self):
        """Close the underlying stream."""
        if not self.stream.closed:
            self.stream.close()


    def flush(self):
        """Flush the underlying stream"""
        if not self.rdonly:
            self.stream.flush()


    def additem(self, pmid, date, features):
        """Add a (pmid,date,features) record to the FeatureStream
        
        @param pmid: PubMed ID (string or integer).
        
        @param date: Either (year,month,day), or YYYMMDD integer date for the
        record.
        
        @param features: Numpy array of feature IDs."""
        if features.dtype != self.ftype:
            raise ValueError("Features dtype must be %s not %s" % 
                             (str(self.ftype),str(features.dtype)))
        if isinstance(date, tuple):
            date = Date2Integer(date)
        #logging.debug("PMID=%d, DATE=%d, ALEN=%d", int(pmid), date, len(features))
        self.stream.write(struct.pack("IIH", int(pmid), date, len(features)))
        features.tofile(self.stream)


    def readitem(self, pos=None):
        """Read a feature stream record from the current location.

        @param pos: Seek here first. Very bad if the position is wrong! 
        
        @return: (PubMed ID, YYYYMMDD, feature vector) as (int,int,array) from
        the stream."""
        if pos is not None:
            self.stream.seek(pos)
        head = self.stream.read(4+4+2)
        if len(head) == 0:
            return None
        pmid, date, alen = struct.unpack("IIH", head)
        #logging.debug("PMID=%d, DATE=%d, ALEN=%d", pmid, date, alen)
        if alen > self.max_bytes:
            raise ValueError("Vector too long (%d) due to bad seek.", alen)
        return (pmid, date, nx.fromfile(self.stream, self.ftype, alen))


    def iteritems(self):
        """Iterate over records obtained via of L{readitem}. This is quite slow
        - ScoreCalculator and FeatureCounter in L{mscanner.fastscores} use C
        code to iterate rapidly over the feature stream."""
        item = self.readitem(0)
        while item is not None:
            yield item
            item = self.readitem()
        self.stream.seek(0,2) # EOF


def encode(numbers):
    """Encode a sorted vector using variable byte encoding of the difference
    between successive items. Each gap is encoded as reverse of the byte list
    one would obtain by writing the least 7 bits with high bit set, followed by
    higher groups of 7 bits until no bits remain.  If C{numbers} is not
    sorted, bad things will happen.
    
    @param numbers: Sorted list of numbers. 
    @return: Numpy uint8  vector with the encoding."""
    result = []
    last = 0
    for number in numbers:
        gap = number - last
        last = number
        # Higher bits are are pushed at the front
        bytes = [0x80 | (gap & 0x7f)]
        gap >>= 7
        while gap > 0:
            bytes.insert(0, gap & 0x7f)
            gap >>= 7            
        result.extend(bytes)
    return nx.array(result, nx.uint8)


def decode(bytestream):
    """Decode a sorted vector whose gaps have been variable-byte encoded. 
    To read a gap: push the least 7 bits of current byte into the least 
    bits of the output, and repeat, stopping when the high bit is set.
    
    @param bytestream: Numpy uint8 vector in variable byte encoding. 
    @yield: Sorted iteration of numbers. """
    gap = 0
    last = 0
    for byte in bytestream:
        gap = (gap << 7) | (byte & 0x7f)
        if byte & 0x80:
            last += gap
            yield last
            gap = 0


class EncodedFeatureStream(PlainFeatureStream):
    """A FeatureStream which compresses the feature vectors using variable byte
    encoding.   
    See http://nlp.stanford.edu/IR-book/html/htmledition/variable-byte-codes-1.html.
    """
    
    def additem(self, pmid, date, features):
        """Add a (pmid,date,features) record to the FeatureStream
        
        @param pmid: PubMed ID (string or integer).
        @param date: (year,month,day), or YYYMMDD integer date.
        @param features: Numpy array of feature IDs."""
        if features.dtype != self.ftype:
            raise ValueError("Features dtype must be %s not %s" % 
                             (str(self.ftype),str(features.dtype)))
        if isinstance(date, tuple):
            date = Date2Integer(date)
        #logging.debug("PMID=%d, DATE=%d, LEN=%d", int(pmid), date, len(features))
        data = encode(features)
        self.stream.write(struct.pack("IIH", int(pmid), date, data.nbytes))
        data.tofile(self.stream)


    def readitem(self, pos=None):
        """Read a feature stream record from the current location.

        @param pos: Seek here first. Very bad if the position is wrong! 
        
        @return: (PubMed ID, YYYYMMDD, feature vector) as (int,int,array) from
        the stream."""
        if pos is not None:
            self.stream.seek(pos)
        head = self.stream.read(4+4+2)
        if len(head) == 0:
            return None
        pmid, date, nbytes = struct.unpack("IIH", head)
        #logging.debug("PMID=%d, DATE=%d, LEN=%d", pmid, date, alen)
        if nbytes > self.max_bytes:
            raise ValueError("Vector too long (%d) due to bad seek.", nbytes)
        data = nx.fromfile(self.stream, nx.uint8, nbytes)
        features = nx.array(list(decode(data)), self.ftype)
        return (pmid, date, features)

    
FeatureStream = EncodedFeatureStream
#FeatureStream = PlainFeatureStream