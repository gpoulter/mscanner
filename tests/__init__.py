#!env python

import unittest

def suite():
    modules = [
        "test_article",
        "test_dbexport",
        "test_dbshelve",
        "test_genedrug",
        "test_medline",
        "test_scoring",
        "test_validation",
        "test_xmlparse"
        ]
    tests = unittest.TestSuite()
    for modname in modules:
        tests.addTest(unittest.findTestCases(__import__(modname)))
    return tests

if __name__ == "__main__":
    unittest.main(defaultTest='suite')
