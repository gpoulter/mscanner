from bsddb import db
from path import path
import tempfile
import unittest

from mscanner.support import dbshelve

class ShelfTests(unittest.TestCase):
    """Test for Shelf class
    
    Tests __get__, len, keys, items, values, iterkeys, iteritems, itervalues, __getitem__, __delitem__
    """

    def setUp(self):
        self.db = dbshelve.open(None)

    def testMethods(self):
        d = self.db
        d["A"] = ("A",2)
        d["B"] = ("B",3)
        self.assertEqual(d["A"],("A",2))
        self.assertEqual(len(d), 2)
        self.assertEqual("B" in d, True)
        self.assertEqual(d.keys(), ["B","A"])
        self.assertEqual(d.items(), [ ("B",("B",3)), ("A",("A",2)) ])
        self.assertEqual(d.values(), [ ("B",3), ("A",2) ])
        self.assertEqual(list(d.iterkeys()), ["B","A"])
        self.assertEqual(list(d.iteritems()), [ ("B",("B",3)), ("A",("A",2)),  ])
        self.assertEqual(list(d.itervalues()), [ ("B",3) , ("A",2), ])
        del d["B"]
        self.assertRaises(KeyError, d.__getitem__, "B")
        self.assertRaises(KeyError, d.__delitem__, "B")

class ShelfTnxTests(ShelfTests):
    """Test for Shelf class transactions
    """

    def setUp(self):
        self.envdir = path(tempfile.mkdtemp(prefix="dbshelve-"))
        self.env = db.DBEnv()
        self.env.open(self.envdir, db.DB_INIT_MPOOL|db.DB_INIT_TXN|db.DB_CREATE)
        self.db = dbshelve.open(self.envdir/'dbshelf.db', db.DB_CREATE|db.DB_AUTO_COMMIT, dbenv=self.env)

    def tearDown(self):
        if self.txn is not None:
            self.txn.abort()
        self.db.close()
        self.env.close()
        self.envdir.rmtree(ignore_errors=True)

    def testMethods(self):
        # Test aborting
        self.txn = self.env.txn_begin()
        self.db.set_txn(self.txn)
        ShelfTests.testMethods(self)
        self.txn.abort()
        self.assertEqual(len(self.db), 0)
        # Test committing
        self.txn = self.env.txn_begin()
        self.db.set_txn(self.txn)
        ShelfTests.testMethods(self)
        self.txn.commit()
        self.assertEqual(len(self.db), 1)
        self.txn = None
        

if __name__ == "__main__":
    unittest.main()
