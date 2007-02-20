import tempfile
from path import path
from random import seed
import unittest
import numpy
from validation import Validator, PerformanceStats
import logging

logging.basicConfig(
    level    = logging.DEBUG,
    datefmt  = "%H:%M:%S",
    format   = "%(asctime)-9s %(levelname)-8s %(message)s",
)

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

    def testPartitionSizes(self):
        """Test that partitioning function for cross-validation"""
        starts, sizes = Validator.partitionSizes(10,5)
        self.assert_((starts == [0,2,4,6,8]).all())
        self.assert_((sizes == [2,2,2,2,2]).all())
        starts, sizes = Validator.partitionSizes(33,5)
        self.assert_((starts == [0,7,14,21,27]).all())
        self.assert_((sizes == [7,7,7,6,6]).all())
        
    def testMoveToFront(self):
        """Tests the swapping function doesn't overwrite"""
        data = numpy.arange(10)
        Validator.moveToFront(data, 3, 3)
        correct = numpy.array([3,4,5,0,1,2,6,7,8,9])
        self.assert_(numpy.all(data == correct))

    def testCrossValidMovement(self):
        """Confirm that data is being shuffled correctly during
        validation by adding these lines to Validator.crossValidate:

            print self.pos[:psize], self.pos[psize:]
            print self.neg[:nsize], self.neg[nsize:]
        """
        val = Validator(
            numfeats = 2,
            featdb = {0:[0], 1:[0,1], 2:[0,1], 3:[0], 4:[0], 5:[0], 6:[0], 7:[0], 8:[0], 9:[0], 10:[0,1]},
            pos = numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
            neg = numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
            nfold = 5,
            pseudocount = 0.1,
            randomise = False)
        pscores, nscores = val.validate()

    def testCrossValid(self):
        """Test that cross-validated scores are correctly calculated"""
        val = Validator(
            numfeats = 3,
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 5:[1,2], 6:[1,2], 7:[0,1,2]},
            pos = numpy.array([0, 1, 2, 3]),
            neg = numpy.array([4, 5, 6, 7]),
            nfold = 4,
            pseudocount = 0.1,
            randomise = False)
        pscores, nscores = val.validate()
        cpscores = numpy.array([-2.39789534,  1.03609192,  1.03609192,  3.43398714])
        cnscores = numpy.array([-3.43398714, -1.03609192, -1.03609192,  2.39789534])
        self.assert_(numpy.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(numpy.allclose(nscores,cnscores,rtol=1e-3))

    def testLeaveOutOne(self):
        """Test of leave-out-one cross validation.  Manually calculate
        scores on the articles to see if they are correct"""
        val = Validator(
            numfeats = 3,
            featdb = {0:[0,1,2], 1:[0,1], 2:[0,1], 3:[0,1], 4:[1,2], 5:[1,2], 6:[1,2], 7:[0,1,2]},
            pos = numpy.array([0, 1, 2, 3]),
            neg = numpy.array([4, 5, 6, 7]),
            nfold = 0,
            pseudocount = 0.1,
            randomise = False)
        pscores, nscores = val.validate()
        cpscores = numpy.array([-2.14126396, 1.30037451, 1.30037451, 1.30037451])
        cnscores = numpy.array([-1.30037451, -1.30037451, -1.30037451,  2.14126396])
        self.assert_(numpy.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(numpy.allclose(nscores,cnscores,rtol=1e-3))

if __name__ == "__main__":
    unittest.main()
