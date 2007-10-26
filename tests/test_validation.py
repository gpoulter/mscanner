"""Test suite for mscanner.validation

                               

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

from mscanner.core.FeatureScores import FeatureScores
from mscanner.core.Validator import LeaveOutValidator, CrossValidator
from mscanner.core.PerformanceStats import PerformanceStats

import logging
logging.basicConfig(level=0)



class PerformanceStatsTests(unittest.TestCase):

    def test_PerformanceStats(self):
        p = PerformanceStats(
            pscores = nx.array([1,2,3,4,4,5,5], dtype=nx.float32),
            nscores = nx.array([0,0,1,1,2,4,5], dtype=nx.float32),
            alpha = 0.5)
        import pprint as pp
        print pp.pformat(p.__dict__)



class ValidatorTests(unittest.TestCase):

    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="valid-"))
        print self.prefix


    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)


    def test_make_partitions(self):
        """Test that partitioning function for cross-validation"""
        starts, sizes = CrossValidator.make_partitions(10,5)
        self.assert_((starts == [0,2,4,6,8]).all())
        self.assert_((sizes == [2,2,2,2,2]).all())
        starts, sizes = CrossValidator.make_partitions(33,5)
        self.assert_((starts == [0,7,14,21,27]).all())
        self.assert_((sizes == [7,7,7,6,6]).all())


    def _make_validator(self, featinfo):
        return CrossValidator(
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 
                      5:[1,2], 6:[1,2], 7:[0,1,2]},
            featinfo = featinfo,
            positives = nx.array([0, 1, 2, 3]),
            negatives = nx.array([4, 5, 6, 7]),
            nfolds = 4,
        )
    

    def _check_scores(self, featinfo, cpscores, cnscores):
        val = self._make_validator(featinfo)
        pscores, nscores = val.validate(randomise=False)
        pp.pprint(pscores)
        pp.pprint(cpscores)
        pp.pprint(nscores)
        pp.pprint(cnscores)
        self.assert_(nx.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(nx.allclose(nscores,cnscores,rtol=1e-3))


    def test1_offsetonly(self):
        print "scores_offsetonly"
        self._check_scores(
            FeatureScores([2,5,7], pseudocount = 0.1, make_scores="scores_offsetonly"),
            nx.array([-2.1016295 ,  1.03609192,  1.03609192,  3.1377213 ]),
            nx.array([-3.1377213 , -1.03609192, -1.03609192,  2.1016295 ]))


    def test2_withabsence(self):
        print "scores_withabsence"
        self._check_scores(
            FeatureScores([2,5,7], pseudocount = 0.1, make_scores="scores_withabsence"),
            nx.array([-2.39789534,  2.20616317,  2.20616317,  4.60405827]),
            nx.array([-4.60405827, -2.20616317, -2.20616317,  2.39789534]))


    def test3_newpseudo(self):
        print "scores_newpseudo"
        self._check_scores(
            FeatureScores([2,5,7], pseudocount = 0.1, make_scores="scores_newpseudo"),
            nx.array([-2.39789534,  1.03609192,  1.03609192,  3.43398714]),
            nx.array([-3.43398714, -1.03609192, -1.03609192,  2.39789534]))


    def test4_oldpseudo(self):
        print "scores_oldpseudo"
        self._check_scores(
            FeatureScores([2,5,7], pseudocount = 0.1, make_scores="scores_oldpseudo"),
            nx.array([-2.39789534,  1.03609192,  1.03609192,  3.43398714]),
            nx.array([-3.43398714, -1.03609192, -1.03609192,  2.39789534]))


    def test5_rubin(self):
        print "scores_rubin"
        self._check_scores(
            FeatureScores([2,5,7], make_scores="scores_rubin"),
            nx.array([-17.32206917,   1.09861231,   1.09861231,  18.420681  ]),
            nx.array([-18.420681  ,  -1.09861231,  -1.09861231,  17.32206917]))


    def test_leaveout_validate(self):
        """Test of leave-out-one cross validation.  Manually calculate
        scores on the articles to see if they are correct"""
        val = LeaveOutValidator(
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 
                      5:[1,2], 6:[1,2], 7:[0,1,2]},
            featinfo = FeatureScores([2,5,7], pseudocount = 0.1),
            positives = nx.array([0, 1, 2, 3]),
            negatives = nx.array([4, 5, 6, 7]),
            nfolds = None,
        )
        pscores, nscores = val.validate()
        cpscores = nx.array([-2.14126396, 1.30037451, 1.30037451, 1.30037451])
        cnscores = nx.array([-1.30037451, -1.30037451, -1.30037451,  2.14126396])
        #pp.pprint(pscores)
        #pp.pprint(cpscores)
        #pp.pprint(nscores)
        #pp.pprint(cnscores)
        self.assert_(nx.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(nx.allclose(nscores,cnscores,rtol=1e-3))



if __name__ == "__main__":
    unittest.main()
