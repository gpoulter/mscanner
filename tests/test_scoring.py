"""Test suite for mscanner.scoring

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import numpy as nx
from path import path
import pprint as pp
import tempfile
import unittest

from mscanner.configuration import rc
from mscanner.medline.FeatureDatabase import FeatureDatabase, FeatureStream
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.FeatureScores import FeatureScores, FeatureCounts
from mscanner import cscore, utils


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


    @utils.usetempfile
    def test_cscore(self, citefname):
        """Tests that the cscore program produces the same output as
        iterScores"""
        featscores = nx.array([0.1, 5.0, 10.0, -5.0, -6.0])
        citations = [(4,[4]), (4,[0,1,2]), (1,[0,2,3]), (2,[0,1]), (3,[1,2,3])]
        offset = 5.0
        # Write citations to disk
        fs = FeatureStream(open(citefname, "w"))
        for pmid, feats in citations:
            fs.write(pmid, nx.array(feats, nx.uint16))
        fs.close()    
        # Calculate scores using Python and cscores
        out_pyscore = list(cscore.pyscore_adaptor(
            citefname, len(citations), featscores, offset, 5, 3))
        out_cscore_pipe = list(cscore.cscore_pipe(
            citefname, len(citations), featscores, offset, 5, 3))
        # Compare Python/cscore for equality
        scores_pipe = nx.array(sorted(score for score,pmid in out_cscore_pipe))
        scores_py = nx.array(sorted(score for score,pmid in out_pyscore))
        self.assert_(nx.allclose(scores_pipe, scores_py))
        # Try to test the ctypes version: cscore2
        print scores_py
        print scores_pipe
        if cscore.score != cscore.cscore_dll: return
        out_cscore_dll = list(cscore.cscore_dll(
            citefname, len(citations), featscores, offset, 5, 3))
        scores_dll = nx.array(sorted(score for score,pmid in out_cscore_dll))
        self.assert_(nx.allclose(scores_dll, scores_py))
        print scores_dll



class FeatureScoresTests(unittest.TestCase):
    """Tests of the L{scoring} module"""
    
    def setUp(self):
        self.pfreqs = nx.array([1,2,0])
        self.nfreqs = nx.array([2,1,0])
        self.pdocs = 2
        self.ndocs = 3
        self.featmap = FeatureMapping()

    def test_pseudocount(s):
        """Constant pseudocount"""
        f = FeatureScores(s.featmap, pseudocount=0.1,
                          make_scores="scores_oldbayes_newpseudo")
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        s.assert_(nx.allclose(
            f.scores, nx.array([-0.35894509,  0.93430924,  0.28768207])))
        
    def test_tfidf(s):
        """TFIDF calculation"""
        f = FeatureScores(s.featmap, pseudocount=0.1)
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        print "TFIDF: %s" % (pp.pformat(f.tfidf),)

    def test_background(s):
        """Background frequency pseudocount"""
        s.featmap.numdocs = 10
        s.featmap.counts = [3,2,1]
        f = FeatureScores(s.featmap, pseudocount=None,
                          make_scores="scores_oldbayes_newpseudo")
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        print "Scores (old): %s" % (pp.pformat(f.scores),)
        s.assert_(nx.allclose(
            f.scores, nx.array([-0.28286278,  0.89381787,  0.28768207])))

    def test_new(s):
        """New score calculation method"""
        s.featmap.numdocs = 10
        s.featmap.counts = [3,2,1]
        f = FeatureScores(s.featmap, pseudocount=None)
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        print "Scores (new): %s" % (pp.pformat(f.scores),)
        print "Offset (new): %s" % str(f.offset)

    def test_cutoff(s):
        """Constant pseudocount and cutoff"""
        f = FeatureScores(s.featmap, pseudocount=0.1,
                          make_scores="scores_oldbayes_oldpseudo",
                          get_postmask="make_rare_positives")
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        #print "Scores (old): %s" % (pp.pformat(f.scores_old),)
        s.assert_(nx.allclose(
            f.scores, nx.array([-0.27193372,  1.02132061,  0.0 ])))

    def test_FeatureCounts(self):
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = FeatureCounts(5, featdb, [1,2,3])
        self.assert_(nx.all(counts == [0,1,2,2,1]))



if __name__ == "__main__":
    unittest.main()
