"""Test suite for mscanner.scoring

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import numpy as nx
from path import path
import tempfile
import unittest

from mscanner import scoring
from mscanner.article import Article
from mscanner.featuremap import FeatureMapping
from mscanner.support.utils import usetempfile

class CScoreTests(unittest.TestCase):
    
    @usetempfile
    def test_CScore(self, citefname):
        """Tests that the cscore program produces the same output as
        iterScores"""
        import numpy as nx
        from mscanner.configuration import rc
        from mscanner.scoring import iterCScores, iterScores
        from mscanner.featuredb import FeatureStream
        import struct
        featscores = nx.array([0.1, 5.0, 10.0, -5.0, -6.0])
        citations = [ (4, [4]), (4, [0,1,2]), (1,[0,2,3]), (2,[0,1]), (3,[1,2,3]) ]
        # Write citations to disk
        fs = FeatureStream(file(citefname, "w"))
        for pmid, feats in citations:
            fs.write(pmid, nx.array(feats, nx.uint16))
        fs.close()    
        # Calculate the scores
        out_c = list(iterCScores(rc.cscore_path, citefname, 
                                  len(citations), featscores, 5, None))
        out_old = list(iterScores(citations, featscores, []))
        # Compare the scores for equality
        scores_c = nx.array(sorted(score for score,pmid in out_c))
        scores_old = nx.array(sorted(score for score,pmid in out_old))
        self.assert_(nx.allclose(scores_c, scores_old))

class ScoringTests(unittest.TestCase):
    """Tests for scoring module functions"""
    
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="scoring-"))
        print self.prefix

    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)
        
    def test_FeatureInfo(self):
        pfreqs = nx.array([1,2,0])
        nfreqs = nx.array([2,1,0])
        pdocs = 2
        ndocs = 3
        featmap = FeatureMapping()
        # With constant pseudocount
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=0.1)
        scores = f.getFeatureScores()
        self.assert_(nx.allclose(
            scores, nx.array([-0.35894509,  0.93430924,  0.28768207])))
        # With background pseudocount
        featmap.numdocs = 10
        featmap.counts = [3,2,1]
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=None)
        scores = f.getFeatureScores()
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.28286278,  0.89381787,  0.28768207])))
        # Constant pseudocount and cutoff
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=0.1, 
                                getFrequencies="getProbabilitiesOldBayes",
                                getPostMask="maskRarePositives")
        scores = f.getFeatureScores()
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.27193372,  1.02132061,  0.0 ])))

    def test_countFeatures(self):
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

if __name__ == "__main__":
    unittest.main()
