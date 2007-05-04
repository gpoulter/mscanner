"""Persistent shelf backed by Berkeley DB

open() -- Open/create a Shelf
Shelf -- Class defining a Berkeley DB-backed shelf

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

from bsddb import db
import cPickle, os, unittest
from path import path
from UserDict import DictMixin
import zlib

def open(filename, flags='c', mode=0660, dbenv=None, txn=None, dbname=None, compress=True):
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
    return Shelf(database, txn, compress)

class Shelf(DictMixin):
    """A shelf to hold pickled objects, built upon a bsddb DB object.  It
    automatically pickles/unpickles data objects going to/from the DB.
    """

    def __init__(self, database, txn=None, do_compression=True):
        """Initialise shelf with a db.DB object

        @param do_compression: If True, use zlib to transparently
        compress pickles.
        """
        self.db = database
        self.do_compression = do_compression
        self.set_txn(txn)
        self.compress = lambda x:x
        self.decompress = lambda x:x
        if self.do_compression:
            self.compress = zlib.compress
            self.decompress = zlib.decompress

    def set_txn(self, txn=None):
        """Set the transaction to use for database operations"""
        self.txn = txn

    def close(self):
        """Close the database.  Shelf must not be used after this"""
        self.db.close()

    def __del__(self):
        self.db.close()

    def __len__(self):
        return self.db.stat()["ndata"]

    def __getitem__(self, key):
        v = self.db.get(key, txn=self.txn)
        if v is None: raise KeyError("Key %s not in database" % repr(key))
        return cPickle.loads(self.decompress(v))

    def __setitem__(self, key, value):
        self.db.put(key, self.compress(cPickle.dumps(value, protocol=2)), self.txn)

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
            yield rec[0], cPickle.loads(self.decompress(rec[1]))
            rec = cur.next()
        cur.close()
        
    def __iter__(self):
        cur = self.db.cursor(self.txn)
        rec = cur.first(dlen=0, doff=0)
        while rec is not None:
            yield rec[0]
            rec = cur.next(dlen=0, doff=0)
        cur.close()
