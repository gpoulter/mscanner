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
        self.prefix.rmtree(ignore_errors=True)
        pass

    def test(self):
        val = Validator(
            featmap = {1:("A","mesh"), 2:("B","year"), 3:("C","mesh"), 4:("D","year"), 5:("E","mesh"), 6:("F","mesh")},
            featdb = {1:[1,2,3], 2:[2,3], 3:[1,3], 4:[4,6], 5:[4], 6:[3,4,5], 7:[4,5]},
            posids = set([1, 2, 3]),
            negids = set([4, 5, 6, 7]),
            nfold = 2,
            pseudocount = 0.1,
            daniel = False,
            genedrug_articles = None,
            )
        seed(0)
        self.assertEqual(
            val.partition(set([1,2,3,4,5,6,7,8,9,10]), 3) ,
            [set([9,2,10,6]), set([8,1,7]), set([3,4,5])] )
        pscores, nscores = val.validate()
        val.report(
            pscores,
            nscores,
            self.prefix,
            path("../lib/templates/style.css"))

if __name__ == "__main__":
    unittest.main()
