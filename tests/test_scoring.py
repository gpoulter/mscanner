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

from mscanner.configuration import rc
from mscanner import cscore, featuredb, featuremap, scoring
from mscanner.support.utils import usetempfile


class CScoreModuleTests(unittest.TestCase):
    """Tests of the L{cscore} package"""
    
    def test_ctypes(self):
        if cscore.score != cscore.cscore_dll: return
        from ctypes import cdll, c_int, byref
        lib = cdll.LoadLibrary(rc.cscore_dll)
        output = c_int()
        lib.double_int(2, byref(output))
        self.assertEqual(output.value, 4)
        lib.double_array.argtypes = [ c_int,
            nx.ctypeslib.ndpointer(dtype=nx.int32, ndim=1, flags='CONTIGUOUS') ]
        a = nx.array([1,2,3,4])
        b = nx.array([2,4,6,8])
        lib.double_array(len(a), a)
        self.assert_(nx.allclose(a, b))


    @usetempfile
    def test_cscore(self, citefname):
        """Tests that the cscore program produces the same output as
        iterScores"""
        featscores = nx.array([0.1, 5.0, 10.0, -5.0, -6.0])
        citations = [ (4, [4]), (4, [0,1,2]), (1,[0,2,3]), (2,[0,1]), (3,[1,2,3]) ]
        # Write citations to disk
        fs = featuredb.FeatureStream(open(citefname, "w"))
        for pmid, feats in citations:
            fs.write(pmid, nx.array(feats, nx.uint16))
        fs.close()    
        # Calculate scores using Python and cscores
        out_pyscore = list(cscore.pyscore_adaptor(citefname, len(citations), featscores, 5, 3))
        out_cscore_pipe = list(cscore.cscore_pipe(citefname, len(citations), featscores, 5, 3))
        # Compare Python/cscore for equality
        scores_pipe = nx.array(sorted(score for score,pmid in out_cscore_pipe))
        scores_py = nx.array(sorted(score for score,pmid in out_pyscore))
        self.assert_(nx.allclose(scores_pipe, scores_py))
        # Try to test the ctypes version: cscore2
        print scores_py
        print scores_pipe
        if cscore.score != cscore.cscore_dll: return
        out_cscore_dll = list(cscore.cscore_dll(
            citefname, len(citations), featscores, 5, 3))
        scores_dll = nx.array(sorted(score for score,pmid in out_cscore_dll))
        self.assert_(nx.allclose(scores_dll, scores_py))
        print scores_dll



class ScoringModuleTests(unittest.TestCase):
    """Tests of the L{scoring} module@usetempfile"""
    
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
        featmap = featuremap.FeatureMapping()
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


    def test_count_features(self):
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = scoring.count_features(5, featdb, [1,2,3])
        self.assert_(nx.all(counts == [0,1,2,2,1]))



if __name__ == "__main__":
    unittest.main()
