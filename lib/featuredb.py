"""Persistent assocation of items to list of integer feautres

@author: Graham Poulter
                                 
"""

import struct
from bsddb import db
from array import array

class FeatureDatabase:
    """Persistent mapping from integer key to array objects of numerical values"""
    
    def __init__(self, filename=None, flags='c', mode=0660, dbenv=None, txn=None, dbname=None, value_type='i'):
        """Initialise database

        @param filename: Database file
        @param flags: Opening flags (r,rw,w,c,n)
        @param mode: Numeric file permissions
        @param dbenv: Optional database environment
        @param txn: Optional database transaction
        @param value_type: Typecode for packing/unpacking of value structs, typically 'H' or 'i'
        """
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
        self.value_type = value_type
        self.db.open(filename, dbname, db.DB_HASH, flags, mode, txn=txn)

    def __del__(self):
        if hasattr(self,"db"):
            self.close()

    def close(self):
        """Close the database.  Do not use this object after doing so"""
        self.db.close()
        delattr(self, "db")

    def getitem(self, key, txn=None):
        """Return an array object of values for a given key"""
        values_packed = self.db.get(struct.pack("i",int(key)), txn=txn)
        if values_packed is None:
            raise KeyError("Record %d not found in feature database" % key)
        return array(self.value_type, values_packed)

    def setitem(self, key, values, txn=None):
        """Associate integer key with an array object of values"""
        if not isinstance(values, array):
            values = array(self.value_type, values)
        if values.typecode != self.value_type:
            raise ValueError("Data value type mismatch: %s vs %s" % (values.typecode, self.value_type))
        self.db.put(struct.pack("i",int(key)), values.tostring(), txn=txn)
        
    def delitem(self, key, txn=None):
        """Delete a given key from the database"""
        self.db.delete(struct.pack("i",int(key)), txn=txn)

    # Bunch of dictionary methods

    def __getitem__(self, key):
        return self.getitem(key)

    def __setitem__(self, key, values):
        self.setitem(key, values)

    def __len__(self):
        return self.db.stat()["ndata"]

    def __contains__(self, key):
        return self.db.has_key(struct.pack("i",int(key)))

    def keys(self):
        return [ k for k in self ]

    def __iter__(self):
        cur = self.db.cursor()
        rec = cur.first(dlen=0, doff=0)
        while rec is not None:
            yield struct.unpack("i",rec[0])[0]
            rec = cur.next(dlen=0, doff=0)
        cur.close()

    def iteritems(self):
        cur = self.db.cursor()
        rec = cur.first()
        while rec is not None:
            yield struct.unpack("i",rec[0])[0], array(self.value_type,rec[1])
            rec = cur.next()
        cur.close()

