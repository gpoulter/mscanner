"""Provides a mapping between feature strings and their integer IDs"""

from __future__ import with_statement
from contextlib import closing
from itertools import izip
import logging
import numpy as nx
try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3

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


class FeatureMapping:
    """Persistent mapping between string features and feature IDs

    Database columns are id (feature ID), type (feature type),
    name (feature string), and count (feature count).  
    
    Feature type is one of "mesh", "qual", "issn", "w" (word) and "a" (author).
    Database is indexed by id and (type,name).
    
    @ivar con: The SQLite database connection.
    
    @ivar filename: Path to SQLite database file.
    
    @ivar grow_features: If True, add new feature strings when encountered. If
    False, ignore unknown feature strings (maintains a static feature space).
    """

    def __init__(self, filename, grow_features=True):
        """Initialise the table of feature counts."""
        self.filename = filename
        self.grow_features = grow_features
        self.filename = filename
        self.con = sqlite3.connect(filename or ":memory:")
        #self.con.execute("PRAGMA synchronous=OFF")
        self.con.execute("PRAGMA cache_size=10000")
        self.con.execute("""CREATE TABLE IF NOT EXISTS fmap (
          id INTEGER PRIMARY KEY,
          type TEXT, name TEXT, count INTEGER,
          UNIQUE(type,name) )""")
        # Make sure that feature 0 exists (create a dummy if necessary)
        if self.con.execute("SELECT id FROM fmap WHERE id=0").fetchone() is None:
            self.con.execute("INSERT OR IGNORE INTO fmap VALUES(0, '', '', 0)")


    def close(self):
        """Close the underlying database"""
        self.con.commit()
        self.con.close()


    def commit(self):
        """Commit pending transactions in the underlying connection"""
        self.con.commit()


    def __len__(self):
        """Return number of features in the table."""
        try:
            return self._length
        except AttributeError:
            if hasattr(self, "_counts"):
                self._length = len(self._counts)
            else:
                logging.debug("FeatureMapping: Query number of features.")
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
            self._counts = nx.zeros(len(self), nx.uint32)
            logging.debug("FeatureMapping: Query occurrence counts vector.")
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
        logging.info("FeatureMapping: Keeping %d, deleting %d out of %d (min %d occurrences).",
                     nkeep, numfeats-nkeep, numfeats, mincount)
        # lookup[oldid] == newid or -1
        lookup = nx.zeros(len(self.counts),nx.int32) - 1
        lookup[keep] = nx.arange(0, nkeep, dtype=nx.int32)
        # Map oldfeatures[keep] -> lookup[keep], delete oldfeatures[~keep]
        oldfeatures = nx.arange(0, numfeats, dtype=nx.int32)
        # Perform the feature deletion and update
        self.con.execute("DELETE FROM fmap WHERE count<?", (mincount,))
        for newid, oldid in izip(lookup[keep],oldfeatures[keep]):
            self.con.execute("UPDATE fmap SET id=? WHERE id=?",
                             (int(newid), int(oldid)))
        self.con.commit()
        self.con.execute("VACUUM")
        delattrs(self, "_counts", "_length")
        return lookup


    def get_feature(self, fid):
        """Retrieve feature (name, type) for a given feature ID. 
        @raise KeyError: if feature ID does not exist."""
        row = self.con.execute(
            "SELECT name, type FROM fmap WHERE id=?", (int(fid),)).fetchone()
        if row is None:
            raise KeyError("Invalid key: %s" % str(key))
        return row


    def add_article(self, featuredict):
        """For each feature, insert it with count 1, or increment count of the
        existing feature. Also returns the feature vector (same result as
        make_vector, which does not alter occurrence counts).
        
        @param featuredict: Dictionary keyed by feature type, where each value
        is a list of feature strings of that type C{{'mesh':['A','B']}}."""
        for ftype, featlist in featuredict.iteritems():
            if len(featlist) > 0:
                c = self.con.execute(
                    "UPDATE fmap SET count=(count+1) WHERE type='"+ftype+"' AND name IN "
                    + self.holders(len(featlist)), featlist)
                if c.rowcount < len(featlist) and self.grow_features:
                    c.executemany("INSERT OR IGNORE INTO fmap VALUES(NULL,'"+ftype+"',?,1)",
                                  ((fname,) for fname in featlist))
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




class MemoryFeatureMapping:
    """Implements the same interface as L{FeatureMapping} but does operations
    in-memory instead of with SQL. It can be substituted for greater speed when
    regenerating or adding articles.
    
    @note: Here, C{counts} is a list attribute, not the array property of
    L{FeatureMapping}. The difference breaks scores_bgfreq and query operations.   
    
    We overwrite the whole SQLite database when commit() is called. Memory is
    much faster than SQLite queries for adding articles. However, the data
    structure takes too much memory when word features are used (MeSH features
    are OK), so we must use L{FeatureMapping} in the task queue.

    @ivar filename: The SQLite database to save to.

    @ivar features: List for looking up (fname, ftype) by feature ID.
    
    @ivar counts: List for looking up occurrence count by feature ID.
    
    @ivar feature_ids: Nested dictionary, for looking up feature ID by
    feature type and feature name.
    """
    
    def __init__(self, filename, grow_features=True):
        self.filename = filename
        self.grow_features = grow_features
        # Simulate insertion of the row (0,"","",0)
        self.features = [("","")]
        self.counts = [0]
        self.feature_ids = {"":{"":0}}
        self.load()

    def close(self):
        pass

    def commit(self):
        if self.filename is None:
            return
        if self.filename.exists():
            self.filename.remove()
        with closing(sqlite3.connect(self.filename or ":memory:")) as con:
            con.execute("""CREATE TABLE fmap (
              id INTEGER PRIMARY KEY,
              type TEXT, name TEXT, count INTEGER,
              UNIQUE(type,name) )""")
            for fid, count in enumerate(self.counts):
                fname, ftype = self.features[fid]
                con.execute("INSERT INTO fmap VALUES(?,?,?,?)", 
                            (fid, ftype, fname, count))
            con.commit()

    def load(self):
        if self.filename is None or not self.filename.exists():
            return
        self.features = []
        self.counts = []
        self.feature_ids = {}
        from contextlib import closing
        with closing(sqlite3.connect(self.filename or ":memory:")) as con:
            for fid, ftype, fname, count in con.execute("SELECT * from fmap"):
                self.features.append((ftype, fname))
                self.counts.append(int(count))
                if ftype not in self.feature_ids:
                    self.feature_ids[ftype] = {}
                self.feature_ids[ftype][fname] = fid

    def __len__(self):
        return len(self.counts)

    def get_feature(self, fid):
        return self.features[fid]

    def type_mask(self, ftypes):
        mask = nx.zeros(len(self), nx.bool)
        if len(ftypes) > 0:
            for ftype in ftypes:
                for fid in self.feature_ids[ftype].itervalues():
                    mask[fid] = True
        return mask

    def add_article(self, featuredict):
        for ftype, featlist in featuredict.iteritems():
            if ftype not in self.feature_ids:
                self.feature_ids[ftype] = {}
            fdict = self.feature_ids[ftype]
            for fname in featlist:
                if fname in fdict:
                    self.counts[fdict[fname]] += 1
                elif self.grow_features:
                    fdict[fname] = len(self.features)
                    self.features.append((fname,ftype))
                    self.counts.append(1)

    def make_vector(self, featuredict):
        vector = []
        for ftype, featlist in featuredict.iteritems():
            for fname in featlist:
                vector.append(self.feature_ids[ftype][fname])
        vector.sort() 
        return vector
