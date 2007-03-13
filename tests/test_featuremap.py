from path import path
import numpy as nx
import tempfile
import unittest

from mscanner.featuremap import FeatureMapping

class TempFileTestCase(unittest.TestCase):
    def setUp(self):
        self.fn = path(tempfile.mktemp())
    def tearDown(self):
        self.fn.remove()

class FeatureMappingTests(TempFileTestCase):
    """Test for FeatureMapping class
    """

    def test(self):
        fm = FeatureMapping(self.fn)
        self.assert_(nx.all(fm.addArticle(Q=["A","B"], T=["A","C"]) == [0,1,2,3]))
        self.assertEqual([fm[i] for i in [0,1,2,3,]], [("A","Q"), ("B","Q"),("A","T"),("C","T")])
        self.assertEqual(fm[1], ("B","Q"))
        self.assertEqual(fm[("C","T")], 3)
        self.assert_(nx.all(fm.counts == [1,1,1,1]))
        fm.dump()
        fm.load()
        fm.dump()
        fm.load()
        self.assertEqual(fm.features, [("A","Q"),("B","Q"),("A","T"),("C","T")])
        self.assertEqual(fm.feature_ids, {"Q":{"A":0,"B":1}, "T":{"A":2,"C":3}})
        self.assert_(nx.all(fm.featureTypeMask("Q") == [1,1,0,0]))

if __name__ == "__main__":
    unittest.main()
