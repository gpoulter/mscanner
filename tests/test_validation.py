#!env python

import tempfile
from path import path
from random import seed
import unittest
from validation import Validator

class ValidatorTest(unittest.TestCase):
    """
    Tests Validator: partition, validate, report
    Implicitly tests Validator: plot*, makeHistogram, articleIsPositive
    """
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="valid-"))
        print self.prefix

    def tearDown(self):
        #self.prefix.rmtree(ignore_errors=True)
        pass

    def test(self):
        val = Validator(
            featmap = [("A","mesh"), ("B","year"), ("C","mesh"), ("D","year"), ("E","mesh"), ("F","mesh")],
            featdb = {1:[0,1,2], 2:[1,2], 3:[0,2], 4:[3,5], 5:[3], 6:[2,3,4], 7:[3,4]},
            posids = [1, 2, 3],
            negids = [4, 5, 6, 7],
            nfold = 2,
            pseudocount = 0.1,
            daniel = False,
            genedrug_articles = None,
            )
        seed(0)
        pscores, nscores = val.validate()
        val.report(
            pscores,
            nscores,
            self.prefix,
            path("../lib/templates/style.css"))

if __name__ == "__main__":
    unittest.main()
