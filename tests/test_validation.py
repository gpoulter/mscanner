"""Test suite for mscanner.validation"""

import logging
import numpy as nx
import pprint as pp
import unittest

from mscanner.configuration import rc
from mscanner.core.FeatureScores import FeatureScores
from mscanner.core.Validator import cross_validate, make_partitions, count_features
from mscanner.core.metrics import PerformanceVectors
from utility import start_logger


class PerformanceVectorsTests(unittest.TestCase):

    def test_PerformanceVectors(self):
        p = PerformanceVectors(
            pscores = nx.array([1,2,3,4,4,5,5], dtype=nx.float32),
            nscores = nx.array([0,0,1,1,2,4,5], dtype=nx.float32),
            alpha = 0.5)
        # Don't know how to test this for correctness - it would take
        # a day to evaluate answers by hand.  Rather, just check that it
        # doesn't give errors, and see if results are consistent on real data.
        logging.debug("PerformanceVectors: %s", pp.pformat(p.__dict__))


class ValidatorTests(unittest.TestCase):


    def test_make_partitions(self):
        """Test the calculation of split-points for partitioning the data"""
        starts, sizes = make_partitions(10,5)
        self.assert_((starts == [0,2,4,6,8]).all())
        self.assert_((sizes == [2,2,2,2,2]).all())
        starts, sizes = make_partitions(33,5)
        self.assert_((starts == [0,7,14,21,27]).all())
        self.assert_((sizes == [7,7,7,6,6]).all())


    def test_count_features(self):
        """Count occurrences of features"""
        features = [[1,2], [2,3], [3,4]]
        counts = count_features(5, features)
        self.assert_(nx.all(counts == [0,1,2,2,1]))


    def _check_scores(self, featinfo, cpscores, cnscores):
        positives = [ [1,2,3], [1,3], [1,3], [1,3] ]
        negatives = [ [1,2], [1,2], [1,2], [1,2,3] ]
        pscores, nscores = cross_validate(featinfo, positives, negatives, 4)
        logging.debug("pscores  output: %s", pp.pformat(pscores))
        logging.debug("pscores correct: %s", pp.pformat(cpscores))
        logging.debug("nscores  output: %s", pp.pformat(nscores))
        logging.debug("nscores correct: %s", pp.pformat(cnscores))
        self.assert_(nx.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(nx.allclose(nscores,cnscores,rtol=1e-3))


    def test_cross_validate(self):
        """Cross validation with scores_laplace_split."""
        rc.mincount = 0
        rc.min_infogain = 0
        rc.type_mask = []
        self._check_scores(
            FeatureScores([0,8,5,5], "scores_laplace_split"),
            nx.array([-0.69314718,  1.79175949,  1.79175949,  2.48490667], nx.float32),
            nx.array([-2.48490667, -1.79175949, -1.79175949,  0.69314718], nx.float32))


if __name__ == "__main__":
    start_logger()
    unittest.main()
