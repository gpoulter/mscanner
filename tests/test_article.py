#!env python

from article import *
from path import path
import tempfile
import unittest

class TempFileTestCase(unittest.TestCase):
    def setUp(self):
        self.fn = path(tempfile.mktemp())
    def tearDown(self):
        self.fn.remove()

class FileTrackerTests(TempFileTestCase):
    """Test for FileTracker class

    Tests add, toprocess, dump, __init__
    """
    def test(self):
        t = FileTracker(self.fn)
        t.add(path("hack/a.xml"))
        t.add(path("cough/b.xml"))
        self.assertEqual(t.toprocess([path("foo/a.xml"), path("blah/c.xml")]), ["blah/c.xml"])
        t.dump()
        del t
        t = FileTracker(self.fn)
        self.assertEqual(t, set(['a.xml', 'b.xml']))
        
class FeatureMappingTests(TempFileTestCase):
    """Test for FeatureMapping class

    Tests: getids, getterms, __get__, dump, load
    """
    def test(self):
        fm = FeatureMapping(self.fn)
        self.assertEqual(fm.getFeatureIds(["A","B"]), [0,1])
        self.assertEqual(fm.getFeatureIds(["A","C"],"T"), [2,3])
        self.assertEqual(fm.getFeatures([0,1,2,3]), [("A",None), ("B",None),("A","T"),("C","T")])
        self.assertEqual(fm[1], ("B",None))
        fm.dump()
        fm.load()
        fm.dump()
        fm.load()
        self.assertEqual(fm.feats, [("A",None),("B",None),("A","T"),("C","T")])
        self.assertEqual(fm.feat2id, {None:{"A":0,"B":1}, "T":{"A":2,"C":3}})
        
class TermCountsTests(TempFileTestCase):
    """Test for TermCounts class

    Tests: add, __get__, dump, load, subtract
    """
    def test(self):
        t = TermCounts(self.fn)
        t.add([1,3])
        t.add([2,3])
        self.assertEqual(t[1], 1)
        self.assertEqual(t[2], 1)
        self.assertEqual(t[3], 2)
        self.assertEqual(t.docs, 2)
        self.assertEqual(t.total, 4)
        t.dump()
        del t
        t = TermCounts(self.fn)
        s = TermCounts()
        s.add([1,3])
        r = t.subtract(s)
        self.assertEqual(r[1], 0)
        self.assertEqual(r[2], 1)
        self.assertEqual(r[3], 1)
        self.assertEqual(r.docs, 1)
        self.assertEqual(r.total, 2)

class ArticleTests(unittest.TestCase):
    """Tests for article module functions

    Tests: chooseRandomLines, readPMIDFile, getArticles
    """
    def test(self):
        pass

if __name__ == "__main__":
    unittest.main()
