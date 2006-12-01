#!env python

from featuredb import *
from numpy import *
import unittest

class FeatureDatabaseTests(unittest.TestCase):
    """
    Tests FeatureDatabase getitem, __iter__, iteritems, in, len
    """
    def test( self ):
        d = FeatureDatabase()
        d.setitem(1, [1,3])
        d.setitem(2, [2,3])
        self.assert_(all(d.getitem(1) == array([1,3])))
        self.assert_(all(d.getitem(2) == array([2,3])))
        self.assertRaises(KeyError, d.getitem, 3)
        self.failUnless(1 in d)
        self.failUnless(2 in d)
        self.failIf(3 in d)
        self.assertEqual(d.keys(), ['2','1'])
        self.assertEqual(list(d.__iter__()), ['2','1'])
        self.assertEqual(len(d), 2)
        d.delitem(2)
        self.failIf(2 in d)
        self.assertRaises(ValueError, d.setitem, 4, array([3.3,4]))

if __name__ == "__main__":
    unittest.main()
