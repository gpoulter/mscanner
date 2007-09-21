"""Miscellaneous utility functions and classes"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


def update(instance, variables, exclude=['self']):
    """Update instance attributes from a dictionary
    
    For example, C{update(self, locals())}
    
    @param instance: Instance to update via setattr()
    @param variables: Dictionary of variables
    @param exclude: Variables to exclude, defaults to ['self']
    """
    for k, v in variables.iteritems():
        if k not in exclude:
            setattr(instance, k, v)


def make_random_subset(k, pool, exclude):
    """Choose a random subset of k articles from pool

    This is suitable when pool is large (say, 16 million items), is allowed to
    come out scrambled, and you need to exclude out certain items from
    selection.
    
    @param k: Number of items to choose from pool
    @param pool: Array of items to choose from (will be scrambled!)
    @param exclude: Set of items that may not be chosen
    @return: A new array of chosen items
    """
    from random import randint
    import numpy as nx
    n = len(pool)
    assert 0 <= k <= n
    for i in xrange(k):
        # Non-selected items are in 0 ... n-i-1
        # Selected items are n-i ... n
        dest = n-i-1
        choice = randint(0, dest) # 0 ... n-i-1 inclusive
        while pool[choice] in exclude:
            choice = randint(0, dest)
        # Move the chosen item to the end, where so it will be part of the
        # selected items in the next iteration. Note: this works using single
        # items - it but would break with slices due to their being views into
        # the vector.
        pool[dest], pool[choice] = pool[choice], pool[dest]
    # Phantom iteration: selected are n-k ... n
    return nx.array(pool[n-k:])


def usetempfile(function):
    """Decorator to call a method with a temporary file
    
    We create a temporary file, pass the file path to the decorated function,
    and finally remove the file. This is useful for unit testing."""
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


class FileTransaction(file):
    """Transaction for Cheetah templates to output direct-to-file.
    
    Cheetah defaults to DummyTransaction which creates a huge list and
    joins them up to create a string.  This is way slower than writing to
    file directly.
    
    Usage::
        with FileTransaction("something.html","wb") as ft:
            Template().respond(ft)
    """
    
    def __init__(self, *args, **kw):
        """Open the file, same parameters as for the builtin"""
        file.__init__(self, *args, **kw)
        self.response = self

    def writeln(self):
        """Write a line of output"""
        self.write(txt)
        self.write('\n')

    def getvalue(self):
        """Not implemented"""
        return None

    def __call__(self):
        return self
