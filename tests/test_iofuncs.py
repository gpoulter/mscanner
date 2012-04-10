"""Test suite for mscanner.utils"""

import unittest
from mscanner.core import iofuncs
from utility import usetempfile

class IOTests(unittest.TestCase):

    @usetempfile
    def test_read_pmids(self, fn):
        allpairs = [(10.0,1), (20.0,2), (30.0,3), (40.0,4), (50.0,5)]
        fn.write_lines(["# comment", "1 10", "2 20 blah", "3 30", "4 40", "5 50"])
        includes = [1,2,3,4]
        excludes = [1]
        pmids, broke, excl  = iofuncs.read_pmids_careful(fn, includes, excludes)
        self.assertEqual(list(pmids), [2,3,4])
        self.assertEqual(list(broke), [5])
        self.assertEqual(list(excl), [1])
        pairs = list(iofuncs.read_scores(fn))
        self.assertEqual(pairs, allpairs)


if __name__ == "__main__":
    unittest.main()
