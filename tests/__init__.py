#!/usr/bin/env python

"""Tests of the MScanner modules"""

import unittest

def suite():
    modules = [
        "test_dbshelve",
        "test_featuredb",
        "test_featuremap",
        "test_medline",
        "test_scoring",
        "test_scorefile",
        "test_storage",
        "test_validation"
        ]
    del modules[4]
    tests = unittest.TestSuite()
    for modname in modules:
        tests.addTest(unittest.findTestCases(__import__(modname)))
    return tests

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
