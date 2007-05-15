import numpy as nx
import unittest
import tempfile

from mscanner.utils import *

def usetempfile(function):
    """Given a function taking arguments of 'self' and a path,
    return a function taking only self, with the wrapped function given
    a temporary file"""
    import tempfile
    from path import path
    def tempfile_wrapper(self):
        try:
            fpath = path(tempfile.mktemp())
            return function(self, fpath)
        finally:
            if fpath.isfile():
                fpath.remove()
    return tempfile_wrapper

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
        
    @usetempfile
    def testFileTracker(self, fn):
        """For utils.FileTracker.(__init__, add, toprocess, dump)"""
        t = FileTracker(fn)
        t.add(path("hack/a.xml"))
        t.add(path("cough/b.xml"))
        self.assertEqual(t.toprocess([path("foo/a.xml"), path("blah/c.xml")]), ["blah/c.xml"])
        t.dump()
        del t
        t = FileTracker(fn)
        self.assertEqual(t, set(['a.xml', 'b.xml']))

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
