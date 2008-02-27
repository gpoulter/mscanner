"""Maps PubMed IDs to feature vectors"""

import numpy as nx
import logging
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


from pysqlite2 import dbapi2 as sqlite3
from FeatureStream import vb_encode, vb_decode

class FeatureVectors:
    """Stores documents in feature vector form in an SQLite database. The
    columns are 'pmid' (PubMed id), 'date' (record date) and 'features'
    (variable-byte-encoded vector of feature IDs representing a document).
    
    Note that we do not support removing documents from the database, because
    the feature stream cannot remove documents either.
    
    @ivar con: The SQLite database connection.
    
    @ivar filename: Path to SQLite database file.
    """

    def __init__(self, filename):
        """Initialise the database"""
        self.filename = filename
        if filename is None:
            self.con = sqlite3.connect(":memory:")
        else:
            self.con = sqlite3.connect(filename)
        self.con.execute("""CREATE TABLE IF NOT EXISTS docs (
          pmid INTEGER PRIMARY KEY, date INTEGER, features BLOB)""")


    def close(self):
        """Close the underlying database"""
        self.con.commit()
        self.con.close()


    def __contains__(self, pmid):
        """Detect whether a vector for a document is already stored"""
        return self.con.execute("SELECT pmid FROM docs WHERE pmid=?",
                                (pmid,)).fetchone() is not None

    def __len__(self):
        """Return number of documents in the database. Caches the result
        because it takes a few seconds on 16 million records."""
        try:
            return self._length
        except AttributeError:
            logging.debug("Querying number of documents in database.")
            self._length = self.con.execute("SELECT count(pmid) FROM docs").fetchone()[0]
            return self._length


    def add_record(self, pmid, date, featurevector):
        """Store a feature vector in the databse.
        @param pmid: PubMed ID of document
        @param date: Record date
        @param featurevector: List/array/iterable of feature IDs
        """
        self.con.execute("INSERT INTO docs VALUES(?,?,?)", (pmid, date, 
                          sqlite3.Binary(vb_encode(featurevector))))
        if hasattr(self, "_length"):
            del self._length


    def update_record(self, pmid, date, featurevector):
        """Update a record's vector and date. Parameters as for L{add_record}."""
        self.con.execute("UPDATE docs SET date=?, features=? WHERE pmid=?", 
                         (date, sqlite3.Binary(vb_encode(featurevector)), pmid))
    
    
    def get_records(self, pmidlist):
        """Iterate over (pmid, date, feature vector) for the specified records,
        ordered by increasing PubMed ID. Feature vectors are numpy arrays.
        @param pmidlist: List of PubMed IDs to retrieve.
        """
        for pmid,date,blob in self.con.execute(
            "SELECT * FROM docs WHERE pmid in ("+\
            ",".join(str(x) for x in pmidlist)+")"):
            yield pmid, date, list(vb_decode(blob))


    def create_index(self):
        """Create index on the date column, for fast filtering by date. Do this
        *after* regenerating or parsing from scratch, because the index slows
        down the adding of records."""
        self.con.execute("CREATE INDEX IF NOT EXISTS dateidx ON docs (date)")


    def iteritems(self, mindate=None, maxdate=None):
        """Iterate over (pmid, date, featurevector) from the database.
        @param mindate, maxdate: Limit the dates to consider.
        """
        self.create_index()
        base = "SELECT * FROM docs "
        if mindate is not None and maxdate is not None:
            cursor = self.con.execute(base+"WHERE date BETWEEN ? AND ?", (mindate,maxdate))
        elif mindate is not None:
            cursor = self.con.execute(base+"WHERE date>=?", (mindate,))
        elif maxdate is not None:
            cursor = self.con.execute(base+"WHERE date<=?", (maxdate,))
        else:
            cursor = self.con.execute(base)
        for pmid, date, blob in cursor:
            vector = list(vb_decode(blob))
            yield pmid, date, list(vb_decode(blob))


    def get_random(self, numrecords, mindate=0, maxdate=99999999):
        """Return cursor iterating over PubMed IDs for numrecords random PubMed
        IDs dates in the specified range."""
        self.create_index()
        if numrecords > len(self):
            numrecords = len(self)
        for pmid, date, blob in self.con.execute(
            "SELECT * FROM docs WHERE date BETWEEN ? AND ? ORDER BY random() LIMIT ?",
            (mindate, maxdate, numrecords)):
            yield pmid, date, list(vb_decode(blob))





