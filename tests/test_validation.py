"""Test suite for mscanner.validation

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import logging
import numpy as nx
from path import path
import pprint as pp
import tempfile
import unittest

from mscanner.configuration import rc
from mscanner.core.FeatureScores import FeatureScores
from mscanner.core.Validator import CrossValidator
from mscanner.core.metrics import PerformanceVectors
from mscanner import tests


class PerformanceVectorsTests(unittest.TestCase):

    def test_PerformanceVectors(self):
        p = PerformanceVectors(
            pscores = nx.array([1,2,3,4,4,5,5], dtype=nx.float32),
            nscores = nx.array([0,0,1,1,2,4,5], dtype=nx.float32),
            alpha = 0.5)
        logging.debug("PerformanceVectors: %s", pp.pformat(p.__dict__))
        
    def test_PerformanceRange(self):
        pass



class ValidatorTests(unittest.TestCase):

    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="valid-"))
        logging.debug("Prefix is: %s", self.prefix)


    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)


    def test_make_partitions(self):
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
        pscores, nscores = val.validate(_randomise=False)
        logging.debug("pscores: %s", pp.pformat(pscores))
        logging.debug("pscores should be: %s", pp.pformat(cpscores))
        logging.debug("nscores: %s", pp.pformat(nscores))
        logging.debug("nscores  should be: %s", pp.pformat(cnscores))
        self.assert_(nx.allclose(pscores,cpscores,rtol=1e-3))
        self.assert_(nx.allclose(nscores,cnscores,rtol=1e-3))


    def test_cross_validate(self):
        """Cross validation with scores_laplace_split."""
        rc.mincount = 0
        rc.min_infogain = 0
        rc.type_mask = []
        self._check_scores(
            FeatureScores([2,5,7], "scores_laplace_split"),
            nx.array([-1.09861231,  2.45673585,  2.45673585,  3.55534816]),
            nx.array([-3.55534816, -2.45673585, -2.45673585,  1.09861231]))


if __name__ == "__main__":
    tests.start_logger()
    unittest.main()
