"""Maps PubMed IDs to feature vectors"""

import numpy as nx
import logging
import struct
from mscanner import delattrs


                                     
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

try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3
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


    def commit(self):
        """Commit pending transactions in the underlying connection"""
        self.con.commit()


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
            if hasattr(self, "_pmids"):
                self._length = len(self._pmids)
            else:
                logging.debug("FeatureVectors: Querying number of documents.")
                self._length = self.con.execute("SELECT count(pmid) FROM docs").fetchone()[0]
            return self._length


    def pmids_array(self):
        """Get array of the PubMed IDs in the database, with caching."""
        try:
            return self._pmids
        except AttributeError:
            self._pmids = nx.zeros(len(self), nx.uint32)
            logging.debug("FeatureVectors: Querying for array PubMed IDs.")
            count = 0
            for pmid, in self.con.execute("SELECT pmid FROM docs"):
                self._pmids[count] = pmid
                count += 1
            return self._pmids


    def add_record(self, pmid, date, featurevector):
        """Store a feature vector in the databse.
        @param pmid: PubMed ID of document
        @param date: Record date
        @param featurevector: List/array/iterable of feature IDs, or a string 
        with the vb_encoded representation of the array.
        """
        if not isinstance(featurevector, str):
            featurevector = vb_encode(featurevector)
        self.con.execute("INSERT INTO docs VALUES(?,?,?)", (pmid, date, 
                          sqlite3.Binary(featurevector)))
        delattrs(self, "_length", "_pmids")


    def update_record(self, pmid, date, featurevector):
        """Update a record's vector and date. Parameters are as for L{add_record}."""
        if not isinstance(featurevector, str):
            featurevector = vb_encode(featurevector)
        self.con.execute("UPDATE docs SET date=?, features=? WHERE pmid=?", 
                         (date, sqlite3.Binary(featurevector), pmid))
    
    
    def get_records(self, pmidlist):
        """Iterate over (pmid, date, feature vector) for the specified records,
        ordered by increasing PubMed ID. Feature vectors are numpy arrays.
        @param pmidlist: List of PubMed IDs to retrieve.
        """
        for pmid,date,blob in self.con.execute(
            "SELECT * FROM docs WHERE pmid in ("+\
            ",".join(str(x) for x in pmidlist)+")"):
            yield pmid, date, list(vb_decode(blob))


    def iteritems(self, decode=True):
        """Iterate over (pmid, date, featurevector) from the database. 
        @param decode: If True, featurevector is a list. If False, 
        featurevector is a vb_encoded buffer object."""
        for pmid, date, blob in self.con.execute("SELECT * FROM docs"):
            if decode: 
                blob = list(vb_decode(blob))
            yield pmid, date, blob




def random_subset(k, pool, exclude):
    """Choose a random list of k items from an array
    
    This is a good algorithm when the pool array is large (say, 16 million
    items), we don't mind if the order of pool gets scrambled, and we have to
    exclude certain items from being selected.
    
    @param k: Number of items to choose from pool
    @param pool: Array of items to choose from (will be scrambled!)
    @param exclude: Set of items that may not be chosen
    @return: A new array of the chosen items
    """
    from random import randint
    import numpy as nx
    n = len(pool)
    assert 0 <= k <= n
    for i in xrange(k):
        # Non-selected items are in 0 ... n-i-1
        # Selected items are n-i ... n
        dest = n-i-1
        choice = randint(0, dest) # 0 ... n-i-1 inclusive
        while pool[choice] in exclude:
            choice = randint(0, dest)
        # Move the chosen item to the end, where so it will be part of the
        # selected items in the next iteration. Note: this works using single
        # items - it but would break with slices due to their being views into
        # the vector.
        pool[dest], pool[choice] = pool[choice], pool[dest]
    # Phantom iteration: selected are n-k ... n
    return nx.array(pool[n-k:])

