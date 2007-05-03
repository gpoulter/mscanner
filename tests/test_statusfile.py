import logging as log
import unittest
import tempfile
from path import path

log.basicConfig(level=0)

from mscanner import statusfile

class TestStatusFile(unittest.TestCase):
    
    def setUp(self):
        self.fn = path(tempfile.mktemp())
    
    def tearDown(self):
        if self.fn.isfile():
            self.fn.remove()
    
    def test(self):
        statusfile.update(5)
        statusfile.start(self.fn, 5)
        statusfile.update(3)
        self.assertEqual(self.fn.lines()[1], "3\n")
        statusfile.update(None)
        self.assertEqual(self.fn.lines()[1], "5\n")
        statusfile.read(self.fn)
        self.assertEqual(self.fn.lines()[1], "5\n")
        self.assertEqual(statusfile.progress, 5)
        self.assertRaises(IOError, statusfile.update, 3)
        self.fn.remove()
        statusfile.close()
        statusfile.start(self.fn, 5)
        statusfile.close()
        self.assertEqual(self.fn.exists(), False)

if __name__ == "__main__":
    unittest.main()