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
from mscanner.core.FeatureScores import FeatureScores, FeatureCounts
from mscanner.core.cscore import cscore
from mscanner import tests

class CScoreModuleTests(unittest.TestCase):
    """Tests of the L{cscore} package"""
    
    def test_ctypes(self):
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


    @tests.usetempfile
    def test_cscore(self, docstream):
        """Tests that the cscore program produces the same output as
        iterScores"""
        featscores = nx.array([0.1, 5.0, 10.0, -5.0, -6.0])
        citations = [(1,20010101,[4]), 
                     (2,20020101,[0,1,2]), 
                     (3,20030101,[0,2,3]), 
                     (4,20040101,[0,1]), 
                     (5,20050101,[1,2,3])]
        # Write citations to disk
        docstream = "C:/test"
        fs = FeatureStream(open(docstream, "w"))
        for pmid, date, feats in citations:
            fs.write(pmid, date, nx.array(feats, nx.uint16))
        fs.close()
        scorer = cscore(
            docstream = docstream,
            numdocs = len(citations),
            featscores = featscores,
            offset = 5.0,
            limit = 5,
            threshold = 0.0,
            mindate = 20020101,
            maxdate = 20050101,
            exclude = set([5,8,9]),
            )
        # Calculate scores using Python and cscores
        out_pyscore = list(scorer.pyscore())
        out_pipe = list(scorer.cscore_pipe())
        # Compare Python/cscore for equality
        print "out_py", pp.pformat(out_pyscore)
        print "out_pipe", pp.pformat(out_pipe)
        scores_pipe = nx.array([score for score,pmid in out_pipe])
        scores_py = nx.array([score for score,pmid in out_pyscore])
        self.assert_(nx.allclose(scores_pipe, scores_py))
        # Try to test the ctypes dll version
        out_dll = list(scorer.cscore_dll())
        print "out_dll", pp.pformat(out_dll)
        scores_dll = nx.array([score for score,pmid in out_dll])
        self.assert_(nx.allclose(scores_dll, scores_py))



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
                          make_scores="scores_newpseudo")
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
                          make_scores="scores_newpseudo")
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
                          make_scores="scores_oldpseudo",
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
