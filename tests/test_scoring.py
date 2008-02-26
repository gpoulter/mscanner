"""Test suite for mscanner.scoring

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from __future__ import with_statement
import logging
from contextlib import closing
import numpy as nx
from path import path
import pprint as pp
import tempfile
import unittest

from mscanner.configuration import rc
from mscanner.medline.FeatureDatabase import FeatureDatabase
from mscanner.medline.FeatureStream import FeatureStream
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.core.FeatureScores import FeatureScores, FeatureCounts
from mscanner.fastscores.ScoreCalculator import ScoreCalculator
from mscanner.fastscores.FeatureCounter import FeatureCounter
from mscanner import tests


class CScoreModuleTests(unittest.TestCase):
    """Tests of the L{cscore} package"""
    
    def setUp(self):
        self.citations = [
            (1,20010101,[4]), 
            (2,20020101,[0,1,2]), 
            (3,20030101,[0,2,3]), 
            (4,20040101,[0,1]), 
            (5,20050101,[1,2,3]),
            (6,20050101,[1,140,145])]

    
    """
    def test_ctypes(self):
        try:
            from ctypes import cdll, c_int, byref
        except ImportError:
            return
        lib = cdll.LoadLibrary(ScoreCalculator.score_base+"32.dll")
        output = c_int()
        lib.double_int(2, byref(output))
        self.assertEqual(output.value, 4)
        lib.double_array.argtypes = [ c_int,
            nx.ctypeslib.ndpointer(dtype=nx.int32, ndim=1, flags='CONTIGUOUS') ]
        a = nx.array([1,2,3,4])
        b = a * 2
        lib.double_array(len(a), a)
        self.assert_(nx.allclose(a, b))
    """


    @tests.usetempfile
    def test_FeatureCounter(self, tmpfile):
        """Test the system for fast feature counting"""
        with closing(FeatureStream(tmpfile, rdonly=False)) as fs:
            for pmid, date, feats in self.citations:
                fs.additem(pmid, date, feats)
        fc = FeatureCounter(
            docstream = tmpfile,
            numdocs = len(self.citations),
            numfeats = 150,
            mindate = 20020101,
            maxdate = 20070101,
            exclude = set([4,8,9]))
        p_ndocs, py_counts = fc.py_counts()
        c_ndocs, c_counts = fc.c_counts()
        logging.debug("FeatureCounter.py_counts: %d, %s", p_ndocs, pp.pformat(py_counts))
        logging.debug("FeatureCounter.c_counts: %d, %s", c_ndocs, pp.pformat(c_counts))
        #self.assertEqual(p_ndocs, c_ndocs)
        #self.assert_(nx.allclose(py_counts, c_counts))


    @tests.usetempfile
    def test_ScoreCalculator(self, tmpfile):
        """Consistency test between document score calculators."""
        featscores = nx.array([0.1, 5.0, 10.0, -5.0, -6.0] + [0]*145, nx.float32)
        # Write citations to disk
        with closing(FeatureStream(tmpfile, rdonly=False)) as fs:
            for pmid, date, feats in self.citations:
                fs.additem(pmid, date, feats)
        # Construct the document score calculator
        scorer = ScoreCalculator(
            docstream = tmpfile,
            numdocs = len(self.citations),
            featscores = featscores,
            offset = 5.0,
            limit = 5,
            threshold = 0.0,
            mindate = 20020101,
            maxdate = 20050101,
            exclude = set([5,8,9]))
        # Compare pyscore and cscore_pipe
        out_pyscore = scorer.pyscore()
        out_pipe = scorer.cscore_pipe()
        logging.debug("ScoreCalculator.pyscore: %s", pp.pformat(out_pyscore))
        logging.debug("ScoreCalculator.cscore_pipe: %s", pp.pformat(out_pipe))
        scores_pipe = nx.array([score for score,pmid in out_pipe])
        scores_py = nx.array([score for score,pmid in out_pyscore])
        self.assert_(nx.allclose(scores_pipe, scores_py))
        # Compare pyscore and cscore_dll 
        """
        try: 
            import  ctypes
        except ImportError: 
            return
        out_dll = scorer.cscore_dll()
        logging.debug("out_dll: %s", pp.pformat(out_dll))
        scores_dll = nx.array([score for score,pmid in out_dll])
        self.assert_(nx.allclose(scores_dll, scores_py))
        """


class FeatureScoresTests(unittest.TestCase):
    
    def setUp(s):
        """Put together dummy article database and FeatureMapping"""
        s.articles = [
            # Positive articles
            ["A","B","D"],
            ["A","B"],
            ["A","C"],
            # Negative articles
            ["B","A","E"],
            ["B","A"],
            ["B","C"],
        ]
        s.numdocs = len(s.articles)
        s.pos = range(0,3)
        s.neg = range(3,6)
        # Calculate the necessaries
        s.featmap = FeatureMapping(filename=None)
        for x in s.articles:
            s.featmap.add_article({"Z":x})
        s.featdb = {}
        for i, x in enumerate(s.articles):
            s.featdb[i] = s.featmap.make_vector({"Z":x})
        s.pdocs = len(s.pos)
        s.ndocs = len(s.neg)
        s.pfreqs = FeatureCounts(len(s.featmap), s.featdb, s.pos)
        s.nfreqs = FeatureCounts(len(s.featmap), s.featdb, s.neg)


    def test_TFIDF(s):
        """Calculation of TF-IDF."""
        f = FeatureScores(s.featmap, "scores_laplace_split")
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        #tfidfs = nx.array([])
        #s.assert_(nx.allclose(f.tfidf, tfidfs))
        logging.debug("FeatureScores.tfidf: %s", pp.pformat(f.tfidf))


    def test_InformationGain(s):
        """Calculation of Information Gain."""
        rc.min_infogain = 0.01
        f = FeatureScores(s.featmap, "scores_laplace_split")
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        selected = nx.ones(len(s.featmap), nx.bool)
        selected[0] = False
        IG = f.infogain(selected)
        logging.debug("Information Gain: %s", pp.pformat(IG))
        # Note: first feature (with IG=0) was de-selected
        s.assert_(nx.allclose(IG, 
        [ 0.02592725,  0.02592725,  0.02592725,  0.        ,  0.02592725]))


    def test_FeatureStats(s):
        """Feature statistics calculation."""
        rc.mincount = 2
        f = FeatureScores(s.featmap, "scores_laplace_split")
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        t = f.stats
        s.assertEqual(t.pos_distinct, 3)
        s.assertEqual(t.pos_occurrences, 6)
        logging.debug("Feature Statistics: %s", pp.pformat(t.__dict__))
        logging.debug("P counts: %s", pp.pformat(s.pfreqs))
        logging.debug("N counts: %s", pp.pformat(s.nfreqs))
    

    def test_FeatureCounts(self):
        """FeatureCounts for total occurrences of each feature in a collection."""
        featdb = {1:[1,2], 2:[2,3], 3:[3,4]}
        counts = FeatureCounts(5, featdb, [1,2,3])
        self.assert_(nx.all(counts == [0,1,2,2,1]))


    def _scoremethodtester(s, scoremethod, answer):
        """Test a score calculation method"""
        rc.mincount = 0
        rc.min_infogain = 0
        rc.type_mask = []
        logging.info("FeatureScores(%s)", scoremethod)
        f = FeatureScores(s.featmap, scoremethod)
        f.numdocs = s.numdocs # For scores_bgfreq
        f.update(s.pfreqs, s.nfreqs, s.pdocs, s.ndocs)
        #logging.debug("P freqs: %s", pp.pformat(f.pfreqs))
        #logging.debug("N freqs: %s", pp.pformat(f.nfreqs))
        #logging.debug("Y Success: %s", pp.pformat(f.success_scores))
        #logging.debug("Y Failure: %s", pp.pformat(f.failure_scores))
        logging.debug("Selected: %s", pp.pformat(f.features))
        logging.debug("Base %f, Scores: %s", f.base, pp.pformat(f.scores))
        s.assert_(nx.allclose(f.scores, nx.array(answer)))


    def test_scores_bgfreq(self):
        """Score calculation using background counts."""
        self._scoremethodtester("scores_bgfreq", 
            [ 0, 2.24819092, -2.24819092,  2.248191  ,  0.        , -2.248191  ])


    def test_scores_laplace_split(self):
        """Score calculation using split-laplace counts."""
        self._scoremethodtester("scores_laplace_split", 
            [ 0, 1.43508453, -1.43508453,  1.43508453,  0.        , -1.43508453])

    def test_scores_laplace_split(self):
        """Score calculation using split-laplace counts."""
        self._scoremethodtester("scores_laplace_split", 
            [ 0, 1.43508453, -1.43508453,  1.43508453,  0.        , -1.43508453])




if __name__ == "__main__":
    tests.start_logger()
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(FeatureScoreTests)
    #suite = unittest.TestSuite()
    #suite.addTest(FeatureScoresTests('test_InformationGain'))
    #unittest.TextTestRunner(verbosity=2).run(suite)
