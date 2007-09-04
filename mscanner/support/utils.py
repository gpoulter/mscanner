"""Handle Articles, processed files, features, feature occurrence counts

Module contents:
    - L{selfupdate}: Update instance from method args and locals
    - L{update}: Update instance dictionary
    - L{randomSample}: Choose random items from an array
    - L{preserve_cwd}: Decorator to preserve working directory
    - L{usetempfile}: Decorator to call method with a self-destructing temp file
"""

                                     
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


def selfupdate(onlyargs=False, exclude=[]):
    """Call in any method to set instance attributes from local variables.
    
    @param onlyargs: If True, use only named arguments and not locals
    
    @param exclude: Names of other variables to exclude.

    @note: First argument to the caller is taken to be the "self" to update
    
    @note: Equivalent to 'vars(self).update(vars()); del self.self'
    """
    import inspect
    # Get caller frame (must be disposed of!)
    cf = inspect.currentframe().f_back 
    try:
        # Instance is first argument to the caller
        instance = cf.f_locals[cf.f_code.co_varnames[0]]
        # Get names of variables to use
        if onlyargs:
            varnames = cf.f_code.co_varnames[1:cf.f_code.co_argcount]
        else:
            varnames = cf.f_code.co_varnames[1:]
        # Update instance with those names
        for var in varnames:
            if var not in exclude:
                setattr(instance, var, cf.f_locals[var])
    finally:
        del cf


def update(instance, variables, exclude=['self']):
    """Update instance attributes
    
    @param instance: Instance to update via setattr()
    @param variables: Dictionary of variables
    @param exclude: Variables to exclude, defaults to ['self']
    """
    for k, v in variables.iteritems():
        if k not in exclude:
            setattr(instance, k, v)


def randomSample(k, pool, exclude):
    """Choose a random subset of k articles from pool

    This method is only suitable when pool is large (say, 15 million items),
    can be scrambled, and you need to rule out certain items from being
    selected.
    
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
        # Move the chosen item to the end, where so it will be
        # part of the selected items in the next iteration
        pool[dest], pool[choice] = pool[choice], pool[dest]
        #orig_dest = pool[dest]
        #pool[dest] = pool[choice]
        #pool[choice] = orig_dest
    # Phantom iteration: selected are n-k ... n
    return nx.array(pool[n-k:])


def preserve_cwd(f):
    """Decorator to save working directory and restore it afterwards."""
    def cwd_preserver(*a, **kw):
        import os
        try:
            cwd = os.getcwd()
            return f(*a, **kw)
        finally:
            os.chdir(cwd)
    def decorator_update(dest, src):
        dest.__dict__.update(src.__dict__)
        dest.__doc__ = src.__doc__
        dest.__name__ = src.__name__
        ##dest.__undecorated__ = src # help epydoc someday
    decorator_update(cwd_preserver, f)
    return cwd_preserver


def usetempfile(function):
    """Unittest method decorator. Decorates methods of the form funcname(self,
    filename) so that filename is replaced with a temporary file which is
    deleted at the end of the call. """
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
