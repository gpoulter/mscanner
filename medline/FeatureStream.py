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


class FeatureStream:
    """Class for reading/writing a binary stream of Medline records, consisting
    of PubMed ID, record completion date and a vector of Feature IDs for
    features present in the record.  This stream is 
    
    @ivar stream: File object underlying the stream.
    
    @ivar ftype: Numpy integer type for representing features.
    """
    
    
    max_len = 2000
    """Max feature vector bytes to read (detect stream errors)."""


    def __init__(self, fname, mode="ab", ftype=nx.uint16):
        """Initialise the stream.
        @param fname: 
        @param mode: Open file for reading (rb) or writing (ab).
        """
        if mode not in ["ab","rb"]:
            raise ValueError("Invalid file mode %s" % mode)
        self.stream = open(fname, mode)
        self.ftype = ftype


    def close(self):
        """Close the underlying stream."""
        if not self.stream.closed:
            self.stream.close()


    def flush(self):
        """Flush the underlying stream"""
        self.stream.flush()


    def additem(self, pmid, date, features):
        """Add a (pmid,date,features) record to the FeatureStream
        
        @param pmid: PubMed ID (string or integer).
        
        @param date: Either (year,month,day), or YYYMMDD integer date for the
        record.
        
        @param features: Numpy array of feature IDs."""
        if self.stream.mode != "ab":
            raise NotImplementedError("Stream file is not open for append.")
        if features.dtype != self.ftype:
            raise ValueError("Array dtype must be %s not %s" % 
                             (str(self.ftype),str(features.dtype)))
        if isinstance(date, tuple):
            date = Date2Integer(date)
        #logging.debug("PMID=%d, DATE=%d, ALEN=%d", int(pmid), date, len(features))
        self.stream.write(struct.pack("IIH", int(pmid), date, len(features)))
        features.tofile(self.stream)


    def readitem(self, pos=None):
        """Read a (PubMed ID, YYYYMMDD, feature vector) tuple at the current
        location in the stream. The first two are integers, and the last is a
        numpy integer array.         
        @param pos: Seek here first. Bad things will happen if this is
        wrong!"""
        if pos is not None:
            self.stream.seek(pos)
        head = self.stream.read(4+4+2)
        if len(head) == 0: 
            return None
        pmid, date, alen = struct.unpack("IIH", head)
        #logging.debug("PMID=%d, DATE=%d, ALEN=%d", pmid, date, alen)
        if alen > self.max_len:
            raise ValueError("Feature vector too long (%d) due to bad seek.", alen)
        return (pmid, date, nx.fromfile(self.stream, self.ftype, alen))


    def iteritems(self):
        """Iterate over results of L{readitem}. For iterating over the database
        rather use ScoreCalculator and FeatureCounter in
        L{mscanner.fastscores}."""
        item = self.readitem(0)
        while item is not None:
            yield item
            item = self.readitem()
