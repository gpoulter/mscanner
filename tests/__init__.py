#!/usr/bin/env python

"""Tests of the MScanner modules"""

import unittest
from utility import start_logger


def suite():
    modules = [
        "test_article",
        "test_iofuncs",
        "test_medline",
        "test_scoring",
        "test_shelf",
        "test_storage",
        "test_validation"
        ]
    #del modules[4]
    start_logger()
    tests = unittest.TestSuite()
    for modname in modules:
        tests.addTest(unittest.findTestCases(__import__(modname)))
    return tests


if __name__ == "__main__":
    unittest.main(defaultTest='suite')
