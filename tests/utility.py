"""Testing utility functions"""

import logging

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


def start_logger():
    logging.basicConfig(level=logging.DEBUG,format="%(levelname)-8s %(message)s")



