#!/usr/bin/env python

"""Tests of the MScanner modules"""

import unittest


def usetempfile(function):
    """Decorator to call a method with a temporary file
    
    Create a temporary file, pass the file path to the wrapped function and
    finally remove the file afterwards. Meant for wrapping unit testing methods
    which require access to a temporary file."""
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


def suite():
    modules = [
        "test_shelf",
        "test_medline",
        "test_scoring",
        "test_storage",
        "test_iofuncs",
        "test_validation"
        ]
    del modules[4]
    tests = unittest.TestSuite()
    for modname in modules:
        tests.addTest(unittest.findTestCases(__import__(modname)))
    return tests


if __name__ == "__main__":
    unittest.main(defaultTest='suite')
