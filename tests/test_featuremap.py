"""Test suite for mscanner.featuremap

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from path import path
import numpy as nx
import tempfile
import unittest

from mscanner import featuremap

class FeatureMappingTests(unittest.TestCase):
    """For FeatureMapping class"""

    def setUp(self):
        self.fn = path(tempfile.mktemp())

    def tearDown(self):
        if self.fn.isfile():
            self.fn.remove()

    def test_FeatureMapping(self):
        fm = featuremap.FeatureMapping(self.fn)
        self.assert_(nx.all(fm.add_article(Q=["A","B"], T=["A","C"]) == [0,1,2,3]))
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
        self.assert_(nx.all(fm.get_type_mask("Q") == [1,1,0,0]))


if __name__ == "__main__":
    unittest.main()
