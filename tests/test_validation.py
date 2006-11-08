#!env python

from path import path
import unittest
from validation import Validator
from random import seed

class ValidatorTest(unittest.TestCase):
    """
    Tests Validator: partition, validate, report
    Implicitly tests Validator: plot*, makeHistogram, articleIsPositive
    """
    def test(self):
        val = Validator(
            meshdb = { 1:"A", 2:"B", 3:"C", 4:"D", 5:"E", 6:"F" },
            featdb = { 1:[1,2,3], 2:[2,3], 3:[1,3], 4:[4,6], 5:[4], 6:[3,4,5], 7:[4,5] },
            posids = set([ 1, 2, 3 ]),
            negids = set([ 4, 5, 6, 7]),
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
            path("/tmp"),
            path("../lib/templates/style.css"))

if __name__ == "__main__":
    unittest.main()
