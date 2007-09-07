"""Test suite for mscanner.scorefile

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import unittest
from mscanner.scorefile import readPMIDs, writePMIDScores
from mscanner.support.utils import usetempfile

class ScorefileTests(unittest.TestCase):
    
    @usetempfile
    def testReadPMIDs(self, fn):
        fn.write_lines(["# comment", "1 10", "2 20 blah", "3 30", "4 40", "5 50"])
        includes = [1,2,3,4]
        excludes = [1]
        pmids = list(readPMIDs(fn, includes, excludes, withscores=False))
        self.assertEqual(pmids, [2,3,4])
        pairs = list(readPMIDs(fn, includes, excludes, withscores=True))
        self.assertEqual(pairs, [(20.0,2),(30.0,3),(40.0,4)])

if __name__ == "__main__":
    unittest.main()