"""Persistent assocation of items to list of integer feautres

@author: Graham Poulter
                                 
"""

from bsddb import db
import numpy
import logging as log

class FeatureDatabase:
    """Persistent mapping from integer key to array objects of numerical values"""
    
    def __init__(self, filename=None, flags='c', mode=0660, dbenv=None, txn=None, dbname=None, ftype=numpy.uint16):
        """Initialise database

        @param filename: Database file
        @param flags: Opening flags (r,rw,w,c,n)
        @param mode: Numeric file permissions
        @param dbenv: Optional database environment
        @param txn: Optional database transaction
        @param dbname: Logical database name
        @param ftype: Numpy numeric feature type
        """
        self.ftype = ftype
        if isinstance(flags, basestring):
            if flags == 'r':
                flags = db.DB_RDONLY
            elif flags == 'rw':
                flags = 0
            elif flags == 'w':
                flags = db.DB_CREATE
            elif flags == 'c':
                flags = db.DB_CREATE
            elif flags == 'n':
                flags = db.DB_TRUNCATE | db.DB_CREATE
            else:
                raise db.DBError("Flag %s is not in 'r', 'rw', 'w', 'c' or 'n'"  % str(flags))
        self.db = db.DB(dbenv)
        self.db.open(filename, dbname, db.DB_HASH, flags, mode, txn=txn)

    def __del__(self):
        if hasattr(self,"db"):
            self.close()

    def close(self):
        """Close the database.  Do not use this object after doing so"""
        self.db.close()
        delattr(self, "db")

    def getitem(self, key, txn=None):
        """Return an ndarray object of values for a given key"""
        buf = self.db.get(str(key), txn=txn)
        if buf is None:
            raise KeyError("Record %d not found in feature database" % key)
        return numpy.fromstring(buf, self.ftype)

    def setitem(self, key, values, txn=None):
        """Associate integer key with an ndarray object of values

        If given a list instead of an ndarray, convert to array of
        self.ftype, and make sure values are unique.
        """
        if not isinstance(values, numpy.ndarray):
            values = numpy.unique(numpy.array(values, self.ftype))
        if values.dtype != self.ftype:
            raise ValueError("data value type mismatch: tried to place " + str(values.dtype) + " for key " + str(key))
        try:
            self.db.put(str(key), values.tostring(), txn=txn)
        except ValueError, e:
            log.error("featuredb: Failed to place " + str(key) + " : " +  str(values))
            raise
        
    def delitem(self, key, txn=None):
        """Delete a given key from the database"""
        self.db.delete(str(key), txn=txn)

    # Bunch of dictionary methods

    def __getitem__(self, key):
        return self.getitem(key)

    def __setitem__(self, key, values):
        self.setitem(key, values)

    def __len__(self):
        """Fast way to check number of items in database"""
        return self.db.stat()["ndata"]

    def __contains__(self, key):
        """Test for document ID membership.  Converts ID to a string first."""
        return self.db.has_key(str(key))

    def keys(self):
        return [ k for k in self ]

    def __iter__(self):
        cur = self.db.cursor()
        rec = cur.first(dlen=0, doff=0)
        while rec is not None:
            yield rec[0]
            rec = cur.next(dlen=0, doff=0)
        cur.close()

    def iteritems(self):
        cur = self.db.cursor()
        rec = cur.first()
        while rec is not None:
            yield rec[0], numpy.fromstring(rec[1],self.ftype)
            rec = cur.next()
        cur.close()

