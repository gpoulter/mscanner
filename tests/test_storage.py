"""Test suite for mscanner.utils.Storage"""

import unittest
from mscanner.core.Storage import Storage, RCStorage


class StorageModuleTests(unittest.TestCase):

    def test_Storage(self):
        # Storage
        s1 = Storage(a = 1, b = 2)
        s1.a = 5
        self.assertEqual(s1.a, 5)
        # RCStorage
        s2 = RCStorage(a = 1)
        s2.b = lambda: s2.a + 1
        self.assertEqual(s2.a, 1)
        self.assertEqual(s2.b, 2)


if __name__ == "__main__":
    unittest.main()
