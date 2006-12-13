#!env python

from article import *
import numpy
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
        self.assertEqual(fm.getFeatureIds(["A","B"],"Q",True), [0,1])
        self.assertEqual(fm.getFeatureIds(["A","C"],"T",True), [2,3])
        self.assertEqual(fm.getFeatures([0,1,2,3]), [("A","Q"), ("B","Q"),("A","T"),("C","T")])
        self.assertEqual(fm[1], ("B","Q"))
        self.assert_(numpy.all(fm.freqs == [1,1,1,1]))
        fm.dump()
        fm.load()
        fm.dump()
        fm.load()
        self.assertEqual(fm.feats, [("A","Q"),("B","Q"),("A","T"),("C","T")])
        self.assertEqual(fm.feat2id, {"Q":{"A":0,"B":1}, "T":{"A":2,"C":3}})
        
class ArticleTests(unittest.TestCase):
    """Tests for article module functions

    Tests: countFeatures, readPMIDFile, getArticles
    """
    def test(self):
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = countFeatures(5,featdb,[1,2,3])
        self.assert_(numpy.all(counts == [0,1,2,2,1]))

if __name__ == "__main__":
    unittest.main()
