"""A class for rapid iteration over the records in Medline."""

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


class FeatureStream:
    """Class for reading/writing a binary stream of Medline records, consisting
    of PubMed ID, record completion date and a vector of Feature IDs for
    features present in the record.  This stream is 
    
    @ivar stream: File-like object (read/write/close, binary strings).
    
    @ivar ftype: Numpy integer type for representing features.
    """

    def __init__(self, stream, ftype=nx.uint16):
        self.stream = stream
        self.ftype = ftype


    def close(self):
        """Close the underlying stream."""
        self.stream.close()


    def write(self, pmid, date, features):
        """Add a record to the stream
        
        @param pmid: PubMed ID (string or integer).
        
        @param date: Either (year,month,day), or YYYMMDD integer date for the
        record.
        
        @param features: Numpy array of feature IDs."""
        if features.dtype != self.ftype:
            raise ValueError("Array dtype must be %s not %s" % (
                str(self.ftype),str(features.dtype)))
        if isinstance(date, tuple):
            date = Date2Integer(date)
        self.stream.write(struct.pack("IIH", int(pmid), date, len(features)))
        features.tofile(self.stream)


    def __iter__(self):
        """Iterate over tuples of (PubMed ID, YYYYMMDD, features). The first
        two are integers, and the last is a numpy integer array.
        
        @note: Rather use the C programs for ScoreCalculator and FeatureCounter
        in the L{mscanner.fastscores} package."""
        header_len = 4+4+2 # IIH header
        head = self.stream.read(header_len)
        while len(head) == header_len:
            pmid, date, alen = struct.unpack("IIH", head)
            yield (pmid, date, nx.fromfile(self.stream, self.ftype, alen))
            head = self.stream.read(header_len)


def Date2Integer(date):
    """Given (year,month,day), return the integer representation"""
    return date[0]*10000 + date[1]*100 + date[2]


def Integer2Date(intdate):
    """Given a YYYYMMDD integer, return (year,month,day)"""
    return intdate//10000, (intdate%10000)//100, intdate%100