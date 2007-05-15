import numpy as nx
from path import path
import tempfile
import unittest

from mscanner.article import Article
from mscanner.featuremap import FeatureMapping
from mscanner import scoring

class ScoringTests(unittest.TestCase):
    """Tests for scoring module functions"""
    
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="scoring-"))
        print self.prefix

    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)
        
    def test_FeatureInfo(self):
        """Implicitly tests calculateFeatureScores"""
        pfreqs = nx.array([1,2,0])
        nfreqs = nx.array([2,1,0])
        pdocs = 2
        ndocs = 3
        featmap = FeatureMapping()
        # Without masking of unknown features
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=0.1)
        scores = f.getFeatureScores()
        self.assert_(nx.allclose(
            scores, nx.array([-0.27193372,  1.02132061,  0.37469345])))
        # With background-calculated pseudocounts
        featmap.numdocs = 10
        featmap.counts = [3,2,1]
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=None)
        scores = f.getFeatureScores()
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.24512244,  0.95444249,  0.37469344])))
        # With cutoff (rare positive scoring features set to zero)
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=0.1, cutoff=True)
        scores = f.getFeatureScores()
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.27193372,  1.02132061,  0.0 ])))

    def testCountFeatures(self):
        """For utils.countFeatures()"""
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = scoring.countFeatures(5, featdb, [1,2,3])
        self.assert_(nx.all(counts == [0,1,2,2,1]))

    def test_filterDocuments(self):
        docs = { 1:[0,2], 2:[1,2], 3:[0,1] }
        fscores = nx.array([1.0, 10.0, 100.0])
        self.assertEqual(nx.sum(fscores[docs[1]]), 101.0)
        self.assertEqual(scoring.filterDocuments(
            scoring.iterScores(docs.iteritems(), fscores, []), 2, 0 ), 
            [ (110.0,2), (101.0,1) ])
        self.assertEqual(scoring.filterDocuments(
            scoring.iterScores(docs.iteritems(), fscores, []), 10, 102.0 ), 
            [ (110.0,2) ])

    #def test_writeReport(self):
    #    pfreqs = nx.array([1,1,2])
    #    nfreqs = nx.array([1,1,0])
    #    pdocs = 2
    #    ndocs = 1
    #    articles = {
    #        "1111": Article(1111,"T","A",meshterms=set(["A","B"])),
    #        "2222": Article(2222,"T","A",meshterms=set(["A","C"])),
    #        "3333": Article(3333,"T","A",meshterms=set(["B","C"])),
    #        }
    #    scores = [(3333,3.5), (2222,2.0), (1111,1.0)]
    #    featmap = [("A","T"), ("B","T"), ("C","Q")]

if __name__ == "__main__":
    unittest.main()
