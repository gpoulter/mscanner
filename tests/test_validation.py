import logging
import numpy as nx
from path import path
from random import seed
import tempfile
import unittest

from mscanner.scoring import FeatureInfo
from mscanner.validation import Validator, PerformanceStats

logging.basicConfig(level=0)

class PerformanceStatisticTest(unittest.TestCase):
    """Tests for PerformanceStats class"""
    def test(self):

        p = PerformanceStats(
            pscores = nx.array([-2,0,1,2,3,4,5,6,7], dtype=nx.float32),
            nscores = nx.array([-3,-2,-1,0,4,5], dtype=nx.float32),
            alpha = 0.5
            )

class ValidatorTest(unittest.TestCase):
    """
    Tests for Validator class (and implicitly, that the plotting module
    at least functions)
    """
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="valid-"))
        print self.prefix

    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)

    def testPartitionSizes(self):
        """Test that partitioning function for cross-validation"""
        starts, sizes = Validator.partitionSizes(10,5)
        self.assert_((starts == [0,2,4,6,8]).all())
        self.assert_((sizes == [2,2,2,2,2]).all())
        starts, sizes = Validator.partitionSizes(33,5)
        self.assert_((starts == [0,7,14,21,27]).all())
        self.assert_((sizes == [7,7,7,6,6]).all())
        
    def testCrossValid(self):
        """Test that cross-validated scores are correctly calculated"""
        featinfo = FeatureInfo([2,5,7], pseudocount = 0.1)
        val = Validator(
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 5:[1,2], 6:[1,2], 7:[0,1,2]},
            featinfo = featinfo,
            positives = nx.array([0, 1, 2, 3]),
            negatives = nx.array([4, 5, 6, 7]),
            nfolds = 4,
            alpha = 0.5,
        )
        pscores, nscores = val.crossValidate(randomise=False)
        cpscores = nx.array([-2.39789534,  1.03609192,  1.03609192,  3.43398714])
        cnscores = nx.array([-3.43398714, -1.03609192, -1.03609192,  2.39789534])
        self.assert_(nx.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(nx.allclose(nscores,cnscores,rtol=1e-3))

    def testLeaveOutOne(self):
        """Test of leave-out-one cross validation.  Manually calculate
        scores on the articles to see if they are correct"""
        featinfo = FeatureInfo([2,5,7], pseudocount = 0.1)
        val = Validator(
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 5:[1,2], 6:[1,2], 7:[0,1,2]},
            featinfo = featinfo,
            positives = nx.array([0, 1, 2, 3]),
            negatives = nx.array([4, 5, 6, 7]),
            nfolds = None,
            alpha = 0.5,
        )
        pscores, nscores = val.validate()
        cpscores = nx.array([-2.14126396, 1.30037451, 1.30037451, 1.30037451])
        cnscores = nx.array([-1.30037451, -1.30037451, -1.30037451,  2.14126396])
        self.assert_(nx.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(nx.allclose(nscores,cnscores,rtol=1e-3))

if __name__ == "__main__":
    unittest.main()
