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

def have_cscore2():
    """Return true if we have ctypes library and the cscore2 DLL"""
    try:
        import ctypes
    except ImportError:
        return False
    # What about cscore.so on Unix systems?
    return path("cscore2.dll").exists()
    

class CScoreTests(unittest.TestCase):
    
    def test_ctypes(self):
        if not have_cscore2(): return
        from ctypes import cdll, c_int, byref
        cscore = cdll.LoadLibrary(r"cscore2.dll")
        output = c_int()
        cscore.double_int(2, byref(output))
        self.assertEqual(output.value, 4)
        import numpy as nx
        cscore.double_array.argtypes = [ c_int,
            nx.ctypeslib.ndpointer(dtype=nx.int32, ndim=1, flags='CONTIGUOUS') ]
        a = nx.array([1,2,3,4])
        cscore.double_array(len(a), a)


    @usetempfile
    def test_CScore(self, citefname):
        """Tests that the cscore program produces the same output as
        iterScores"""
        import numpy as nx
        from mscanner.configuration import rc
        from mscanner.scoring import iterScores, iterCScores, iterCScores2
        from mscanner.featuredb import FeatureStream
        import struct
        featscores = nx.array([0.1, 5.0, 10.0, -5.0, -6.0])
        citations = [ (4, [4]), (4, [0,1,2]), (1,[0,2,3]), (2,[0,1]), (3,[1,2,3]) ]
        # Write citations to disk
        fs = FeatureStream(open(citefname, "w"))
        for pmid, feats in citations:
            fs.write(pmid, nx.array(feats, nx.uint16))
        fs.close()    
        # Calculate scores using Python and cscores
        out_cscore = list(iterCScores(
            rc.cscore_path, citefname, len(citations), featscores, 5, 3))
        out_old = list(iterScores(citations, featscores, 5))
        # Compare Python/cscore for equality
        scores_cscore = nx.array(sorted(score for score,pmid in out_cscore))
        scores_old = nx.array(sorted(score for score,pmid in out_old))
        print scores_old
        print scores_cscore
        self.assert_(nx.allclose(scores_cscore, scores_old))
        # Try to test the ctypes version: cscore2
        if not have_cscore2(): return
        out_cscore2 = list(iterCScores2(
            citefname, len(citations), featscores, 5, 3))
        scores_cscore2 = nx.array(sorted(score for score,pmid in out_cscore2))
        print scores_cscore2
        self.assert_(nx.allclose(scores_cscore2, scores_old))



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
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.35894509,  0.93430924,  0.28768207])))
        print f.tfidf
        # With background pseudocount
        featmap.numdocs = 10
        featmap.counts = [3,2,1]
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=None)
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.28286278,  0.89381787,  0.28768207])))
        # Constant pseudocount and cutoff
        f = scoring.FeatureInfo(featmap, pfreqs, nfreqs, pdocs, ndocs, 
                                pseudocount=0.1, 
                                frequency_method="getProbabilitiesOldBayes",
                                post_masker="maskRarePositives")
        self.assert_(nx.allclose(
            f.scores, nx.array([-0.27193372,  1.02132061,  0.0 ])))


    def test_countFeatures(self):
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = scoring.countFeatures(5, featdb, [1,2,3])
        self.assert_(nx.all(counts == [0,1,2,2,1]))



if __name__ == "__main__":
    unittest.main()
