"""Provides a mapping between feature strings and their integer IDs"""

from mscanner import delattrs
import logging
import numpy as nx
import sqlite3

                                     
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
        # Make sure a dummy feature occupies the zeroth index
        if self.con.execute("SELECT * FROM fmap LIMIT 1").fetchone() is None:
            self.con.execute("INSERT INTO fmap VALUES(0, '', '', 0)")


    def close(self):
        """Close the underlying database"""
        self.con.commit()
        self.con.close()


    def get_feature(self, fid):
        """Retrieve (feature name, feature type) given feature ID."""
        row = self.con.execute(
            "SELECT name, type FROM fmap WHERE id=?", (fid,)).fetchone()
        if row is None:
            raise KeyError("Invalid key: %s" % str(key))
        return row


    def __len__(self):
        """Return number of features in the table"""
        return self.con.execute("SELECT count(id) FROM fmap").fetchone()[0]


    @staticmethod
    def holders(n):
        """Return a '(?,?,?)' place-holder tuple with n question marks"""
        if n < 1:
            raise ValueError("Must have at least one '?' in the tuple")
        return "(?" + (",?" * (n-1)) + ")"


    def type_mask(self, ftypes):
        """Boolean array to mark features of the specified type.
        
        @param ftypes: List of feature types to exclude, e.g. ["mesh","issn"].
        
        @return: Boolean array indexed by feature ID. True positions correspond
        to features of the types specified in L{ftypes}."""
        mask = nx.zeros(len(self), nx.bool)
        mask[0] = True # Always mask the dummy feature
        if len(ftypes) > 0:
            for fid, in self.con.execute(
                "SELECT id FROM fmap WHERE type IN "+
                self.holders(len(ftypes)), ftypes):
                mask[fid] = True
        return mask
    
    
    def make_vector(self, featuredict):
        """Calculate vector of feature IDs representing an instance, given a
        dictionary with the types and names of the features of the instance
        
        @param featuredict: Dictionary whose keys are feature types, with
        values being a list of feature names of that type, such as
        C{{'mesh':['A','B','C']}}.
     
        @return: Vector of feature IDs.
        """
        vector = []
        for ftype, featlist in featuredict.iteritems():
            if len(featlist) > 0:
                for fid, in self.con.execute(
                    "SELECT id FROM fmap WHERE type='" + ftype + "' AND name IN "
                    + self.holders(len(featlist)), featlist):
                    vector.append(fid)
        # Sort vector prior to compression by EncodedFeatureStream
        vector.sort() 
        return vector


    @property
    def counts(self):
        """Array with the number of occurrences of each feature. Array index is
        the feature ID.
        
        Note that the first element (0) is a dummy feature with count zero,
        because the feature ID is the SQLite ROWID, which starts from 1."""
        try:
            return self._counts
        except AttributeError:
            self._counts = nx.zeros(len(self), nx.uint32)
            for id,count in self.con.execute("SELECT id,count FROM fmap ORDER BY id"):
                self._counts[id] = count
            return self._counts


    def remove_vector(self, featurevector):
        """Remove an article's feature vector from the feature map by
        decrementing the stored count of each feature (also decrements
        L{counts}). Does NOT remove features whose count has dropped to
        zero."""
        self.con.execute("UPDATE fmap SET count=(count-1) WHERE id IN "+
                         self.holders(len(featurevector)), featurevector)
        self._counts[featurevector] -= 1


    def add_article(self, featuredict):
        """Add an article to the feature map. Increments the occurrence
        count for features already represented.
        
        @param featuredict: Dictionary keyed by feature type, where each value
        is a list of feature strings of that type C{{'mesh':['A','B']}}."""
        c = self.con.cursor()
        for ftype, featurelist in featuredict.iteritems():
            for fname in featurelist:
                row = c.execute("SELECT id FROM fmap WHERE type=? AND name=?", 
                                (ftype, fname)).fetchone()
                if row is None:
                    c.execute("INSERT INTO fmap (type,name,count) VALUES(?,?,?)",
                              (ftype, fname, 1))
                else:
                    c.execute("UPDATE fmap SET count=(count+1) WHERE id=?", row)
        c.close()
        # Will need to recalculate count vector if number of features changed
        delattrs(self, "_counts")
