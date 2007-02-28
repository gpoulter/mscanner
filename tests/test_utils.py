import numpy as nx
from path import path
import tempfile
import unittest

from utils import *

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
        
       
class ArticleTests(TempFileTestCase):
    """Tests for article module functions

    Tests: countFeatures, readPMIDFile, getArticles
    """
    def testCountFeatures(self):
        self.fn.touch()
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = countFeatures(5, featdb, [1,2,3])
        self.assert_(nx.all(counts == [0,1,2,2,1]))

    def testReadPMIDFile(self):
        self.fn.write_lines(["# comment", "1 10", "2 20 blah", "3 30", "4 40", "5 50"])
        includes = [1,2,3,4]
        excludes = [1]
        pmids = list(readPMIDs(self.fn, includes, excludes, withscores=False))
        self.assertEqual(pmids, [2,3,4])
        pairs = list(readPMIDs(self.fn, includes, excludes, withscores=True))
        self.assertEqual(pairs, [(2,20.0),(3,30.0),(4,40.0)])

if __name__ == "__main__":
    unittest.main()
