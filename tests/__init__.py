#!/usr/bin/env python

"""Tests of the MScanner modules"""

import unittest

def suite():
    modules = [
        "test_shelf",
        "test_medline",
        "test_scoring",
        "test_storage",
        "test_utils",
        "test_validation"
        ]
    del modules[4]
    tests = unittest.TestSuite()
    for modname in modules:
        tests.addTest(unittest.findTestCases(__import__(modname)))
    return tests

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
