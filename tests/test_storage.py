import unittest
from mscanner.storage import *

class StorageModuleTests(unittest.TestCase):
    
    def testStorage(self):
        """For Storage, RCStorage and RCStorageProxy"""
        # Storage
        s1 = Storage(a = 1, b = 2)
        s1.a = 5
        self.assertEqual(s1.a, 5)
        # RCStorage
        s2 = RCStorage(a = 1)
        s2.b = lambda: s2.a + 1
        self.assertEqual(s2.a, 1)
        self.assertEqual(s2.b, 2)
        # RCStorageProxy
        s3 = RCStorageProxy(s2)
        self.assertEqual(s3.b, 2)
        s3.a = 3
        self.assertEqual(s2.a, 1)
        self.assertEqual(s3.b, 2)
        self.assertEqual(s3.a, 3)
        # Swap s2 <--> s3
        q = s2
        s2 = s3
        s3 = q
        self.assertEqual(s2.a, 3)
        # since s2.a == 3
        self.assertEqual(s2.b, 4)

if __name__ == "__main__":
    unittest.main()