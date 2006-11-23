"""Persistent shelf backed by Berkeley DB

@author: Graham Poulter
                                   

open() -- Open/create a Shelf
Shelf -- Class defining a Berkeley DB-backed shelf

"""

import cPickle, os, unittest
from path import path
from bsddb import db
from UserDict import DictMixin

def open(filename, flags='c', mode=0660, dbenv=None, dbname=None, txn=None):
    """Open a shelf with Berkeley DB backend

    flag is one of 'r','rw','w','c','n'.  Optionally specify flags
    such as db.DB_CREATE, db.DB_TRUNCATE, db.DB_RDONLY,
    db.DB_AUTO_COMMIT.  dbenv provides a db.DBEnv environment, and
    dbname selects a sub-database from the file.
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
    database = db.DB(dbenv)
    database.open(filename, dbname, db.DB_HASH, flags, mode, txn=txn)
    return Shelf(database, txn)

class Shelf(DictMixin):
    """A shelf to hold pickled objects, built upon a bsddb DB object.  It
    automatically pickles/unpickles data objects going to/from the DB.
    """

    def __init__(self, database, txn=None):
        """Initialise shelf with a db.DB object"""
        self.db = database
        self.set_txn(txn)

    def set_txn(self, txn=None):
        """Set the transaction to use for database operations"""
        self.txn = txn

    def close(self):
        """Close the database.  Shelf must not be used after this"""
        self.db.close()

    def __del__(self):
        self.db.close()

    def __len__(self):
        return len(self.db)

    def __getitem__(self, key):
        v = self.db.get(key, txn=self.txn)
        if v is None: raise KeyError("Key %s not in database" % str(key))
        return cPickle.loads(v)

    def __setitem__(self, key, value):
        self.db.put(key, cPickle.dumps(value, protocol=2), self.txn)

    def __delitem__(self, key):
        self.db.delete(key, self.txn)

    def keys(self):
        return self.db.keys(self.txn)

    def items(self):
        return list(self.iteritems())

    def values(self):
        return [ v for k,v in self.iteritems() ]

    def __contains__(self, key):
        return self.db.has_key(key, self.txn)

    def iteritems(self):
        cur = self.db.cursor(self.txn)
        rec = cur.first()
        while rec is not None:
            yield rec[0], cPickle.loads(rec[1])
            rec = cur.next()
        cur.close()
        
    def __iter__(self):
        cur = self.db.cursor(self.txn)
        rec = cur.first(dlen=0, doff=0)
        while rec is not None:
            yield rec[0]
            rec = cur.next(dlen=0, doff=0)
        cur.close()
