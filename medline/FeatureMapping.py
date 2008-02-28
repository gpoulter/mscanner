"""Provides a mapping between feature strings and their integer IDs"""

from mscanner import delattrs
from itertools import izip
import logging
import numpy as nx
from pysqlite2 import dbapi2 as sqlite3

                                     
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


class FeatureMapping:
    """Persistent mapping between string features and feature IDs

    Database columns are id (feature ID), type (feature type),
    name (feature string), and count (feature count).  
    
    Feature type is one of "mesh", "qual", "issn", "w" (word) and "a" (author).
    Database is indexed by id and (type,name).
    
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
        self.con.execute("""CREATE TABLE IF NOT EXISTS fmap (
          id INTEGER PRIMARY KEY,
          type TEXT, name TEXT, count INTEGER,
          UNIQUE(type,name) )""")
        # Make sure that feature 0 exists (create a dummy if necessary)
        self.con.execute("INSERT OR IGNORE INTO fmap VALUES(0, '', '', 0)")


    def close(self):
        """Close the underlying database"""
        self.con.commit()
        self.con.close()


    def __len__(self):
        """Return number of features in the table."""
        try:
            return self._length
        except AttributeError:
            logging.debug("Querying length of FeatureMapping.")
            self._length = self.con.execute("SELECT count(id) FROM fmap").fetchone()[0]
            return self._length


    @property
    def counts(self):
        """Array with the number of occurrences of each feature, indexed by
        feature ID. First element (index 0) is a dummy feature with count zero,
        because the feature ID/index is the SQLite ROWID, which by default
        starts from 1."""
        try:
            return self._counts
        except AttributeError:
            logging.debug("Querying occurrence counts from FeatureMapping.")
            self._counts = nx.zeros(len(self), nx.uint32)
            for id,count in self.con.execute("SELECT id,count FROM fmap ORDER BY id"):
                self._counts[id] = count
            return self._counts


    @staticmethod
    def holders(n):
        """Return a '(?,?,?)' place-holder tuple with n question marks"""
        if n < 1:
            raise ValueError("Must have at least one '?' in the tuple")
        return "(?" + (",?" * (n-1)) + ")"


    def type_mask(self, ftypes):
        """Get boolean array to mark features of the specified type.
        @param ftypes: List of feature types, such as ["mesh","issn"].
        @return: Boolean array indexed by feature ID where True positions
        correspond to features of the types specified in L{ftypes}."""
        mask = nx.zeros(len(self), nx.bool)
        if len(ftypes) > 0:
            for fid, in self.con.execute(
                "SELECT id FROM fmap WHERE type IN "+
                self.holders(len(ftypes)), ftypes):
                mask[fid] = True
        return mask
    
    
    def remove_vector(self, featurevector):
        """Remove an instance from the background count by decrementing the
        occurrence count of each of its features. Does NOT delete features
        whose count has dropped to zero."""
        self.con.execute("UPDATE fmap SET count=(count-1) WHERE id IN "+
                         self.holders(len(featurevector)), featurevector)
        self._counts[featurevector] -= 1


    def vacuum(self, mincount):
        """Delete features with fewer than the specified number of occurrences.
        @return: Array mapping old feature IDs to new, or -1 for deleted features."""
        keep = self.counts >= mincount
        nkeep, numfeats = sum(keep), len(keep)
        logging.info("Keeping %d features (at least %d occurrences). Deleting %d features out of %d.",
                     nkeep, mincount, numfeats-nkeep, numfeats)
        # lookup[oldid] == newid or -1
        lookup = nx.zeros(len(self.counts),nx.int32) - 1
        lookup[keep] = nx.arange(0, nkeep, dtype=nx.int32)
        # Map oldfeatures[keep] -> lookup[keep], delete oldfeatures[~keep]
        oldfeatures = nx.arange(0, numfeats, dtype=nx.int32)
        # Perform the feature deletion and update
        self.con.execute("DELETE FROM fmap WHERE count<?", (mincount,))
        self.con.executemany("UPDATE fmap SET id=? WHERE id=?", 
                             izip(lookup[keep],oldfeatures[keep]))
        delattrs(self, "_counts", "_length")
        return lookup


    def get_feature(self, fid):
        """Retrieve feature (name, type) for a given feature ID. 
        @raise KeyError: if feature ID does not exist."""
        row = self.con.execute(
            "SELECT name, type FROM fmap WHERE id=?", (fid,)).fetchone()
        if row is None:
            raise KeyError("Invalid key: %s" % str(key))
        return row


    def add_article(self, featuredict):
        """For each feature, insert it with count 1, or increment count of the
        existing feature. Also returns the feature vector (same result as
        make_vector, which does not alter occurrence counts).
        
        @param featuredict: Dictionary keyed by feature type, where each value
        is a list of feature strings of that type C{{'mesh':['A','B']}}."""
        c = self.con.cursor()
        for ftype, featurelist in featuredict.iteritems():
            for fname in featurelist:
                try:
                    c.execute("INSERT INTO fmap VALUES(NULL,?,?,?)", (ftype, fname, 1))
                except sqlite3.IntegrityError:
                    c.execute("UPDATE fmap SET count=(count+1) WHERE type=? AND name=?", (ftype,fname))
        c.close()
        # Feature counts have changed
        delattrs(self, "_counts", "_length")


    def make_vector(self, featuredict):
        """Get array of feature IDs representing an instance, given a
        dictionary with the types and names of the features of the instance
        @param featuredict: Dictionary keyed by feature type, where each value
        is a list of feature strings of that type C{{'mesh':['A','B']}}
        @return: Array of feature IDs. 
        """
        vector = [] # Feature vector
        for ftype, featlist in featuredict.iteritems():
            if len(featlist) > 0:
                for fid, in self.con.execute(
                    "SELECT id FROM fmap WHERE type='" + ftype + "' AND name IN "
                    + self.holders(len(featlist)), featlist):
                    vector.append(fid)
        # Sorted vector needed for variable byte encoding
        vector.sort() 
        return vector
