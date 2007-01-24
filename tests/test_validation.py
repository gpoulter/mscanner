#!env python

import tempfile
from path import path
from random import seed
import unittest
import numpy
from validation import Validator, PerformanceStats

class PerformanceStatisticTest(unittest.TestCase):

    def test(self):
        p = PerformanceStats(
            pscores = numpy.array([-2,0,1,2,3,4,5,6,7], dtype=numpy.float32),
            nscores = numpy.array([-3,-2,-1,0,4,5], dtype=numpy.float32),
            alpha = 0.5
            )

class ValidatorTest(unittest.TestCase):
    """
    Tests Validator: partition, validate, report
    Implicitly tests Validator: plot*, makeHistogram, articleIsPositive
    """
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="valid-"))
        print self.prefix

    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)
        pass

    def testPartitionSizes(self):
        """Test that partitioning function for cross-validation"""
        starts, sizes = Validator.partitionSizes(10,5)
        self.assert_((starts == [0,2,4,6,8]).all())
        self.assert_((sizes == [2,2,2,2,2]).all())
        starts, sizes = Validator.partitionSizes(33,5)
        self.assert_((starts == [0,7,14,21,27]).all())
        self.assert_((sizes == [7,7,7,6,6]).all())

    def testCrossValidMovement(self):
        """Confirm that data is being shuffled correctly during
        validation by adding these lines to Validator.crossValidate:

            print self.pos[:psize], self.pos[psize:]
            print self.neg[:nsize], self.neg[nsize:]
        """
        return
        val = Validator(
            featmap = [("A","X")],
            featdb = {0:[0], 1:[0], 2:[0], 3:[0], 4:[0], 5:[0], 6:[0], 7:[0], 8:[0], 9:[0]},
            pos = numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            neg = numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            nfold = 5,
            pseudocount = 0.1,
            norandom = True)
        pscores, nscores = val.validate()

    def testCrossValid(self):
        """Test that cross-validated scores are correctly calculated"""
        val = Validator(
            featmap = [("A","X"), ("B","Y"), ("C","X")],
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 5:[1,2], 6:[1,2], 7:[0,1,2]},
            pos = numpy.array([0, 1, 2, 3]),
            neg = numpy.array([4, 5, 6, 7]),
            nfold = 4,
            pseudocount = 0.1,
            norandom = True)
        pscores, nscores = val.validate()
        print pscores
        print nscores
        #val.report(pscores, nscores, self.prefix, path("../lib/templates/style.css"))

    def testLeaveOutOne(self):
        """Test of leave-out-one cross validation.  Manually calculate
        scores on the articles to see if they are correct"""
        val = Validator(
            featmap = [("A","X"), ("B","Y"), ("C","X")],
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 5:[1,2], 6:[1,2], 7:[0,1,2]},
            pos = numpy.array([0, 1, 2, 3]),
            neg = numpy.array([4, 5, 6, 7]),
            nfold = 0,
            pseudocount = 0.1,
            norandom = True)
        pscores, nscores = val.validate()
        print pscores
        print nscores

if __name__ == "__main__":
    unittest.main()
