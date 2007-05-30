from numpy import *
from path import path
import tempfile
import unittest

from mscanner.featuredb import FeatureDatabase, FeatureStream

class FeatureDatabaseTests(unittest.TestCase):
    """
    Tests FeatureDatabase getitem, __iter__, iteritems, in, len
    """
    def test(self):
        d = FeatureDatabase()
        d.setitem(1, array([1,3],uint16)) #eliminate duplicate features
        d.setitem(2, array([2,3],uint16))
        self.assert_(all(d.getitem(1) == [1,3]))
        self.assert_(all(d.getitem(2) == [2,3]))
        self.assertRaises(KeyError, d.getitem, 3)
        self.failUnless(1 in d)
        self.failUnless(2 in d)
        self.failIf(3 in d)
        self.assertEqual(d.keys(), ['2','1'])
        self.assertEqual(list(d.__iter__()), ['2','1'])
        self.assertEqual(len(d), 2)
        d.delitem(2)
        self.failIf(2 in d)
        self.assertRaises(ValueError, d.setitem, 4, array([3.3,4]))

class FeatureStreamTests(unittest.TestCase):
    """
    Test FeatureStream
    """
    def setUp(self):
        self.fn = path(tempfile.mktemp())
        
    def tearDown(self):
        if self.fn.isfile():
            self.fn.remove()
        
    def test(self):
        f = file(self.fn, "ab")
        fs = FeatureStream(f)
        pmids = (12,34,56)
        feats = [array([1,2,3,4],uint16), array([5,6,7,8],uint16), array([],uint16)]
        for pmid, feat in zip(pmids,feats):
            fs.write(pmid, feat)
        f.close()
        f = file(self.fn, "rb")
        fs = FeatureStream(f)
        rpmids, rfeats = zip(*[x for x in fs])
        self.assertEqual(pmids, rpmids)
        for a, ra in zip(feats, rfeats):
            self.assert_(all(a == ra))
        f.close()

if __name__ == "__main__":
    unittest.main()
