import numpy as nx
import os
import unittest
import tempfile

from mscanner.support.utils import preserve_cwd, selfupdate

class UtilsTests(unittest.TestCase):
    
    def testSelfUpdate(self):
        """For utils.selfupdate()"""
        class X:
            a = 1
            def __init__(s, a, b):
                c = 3
                selfupdate()
        x = X(2,4)
        self.assertEqual(x.a,2)
        self.assertEqual(x.b,4)
        self.assertEqual(x.c,3)
    
    def testDirectoryPreserver(self):
        """For utils.preserve_cwd()"""
        origwd = os.getcwd()
        tmpdir = tempfile.gettempdir()
        @preserve_cwd
        def dirchange1():
            os.chdir(tmpdir)
            assert os.getcwd() == tmpdir
        dirchange1()
        assert os.getcwd() == origwd
        

if __name__ == "__main__":
    unittest.main()
