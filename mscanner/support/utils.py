"""Handle Articles, processed files, features, feature occurrence counts

preserve_cwd -- Decorator to preserve working directory
selfupdate() -- Update instance from method args and locals
update() -- Update instance dictionary
usetempfile -- Decorator to call method with a self-destructing temp file

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

from mscanner.support import dbshelve

def preserve_cwd(f):
    """Decorator which saves working directory before callijng the function,
    and restores it afterwards."""
    def cwd_preserver(*a, **kw):
        import os
        try:
            cwd = os.getcwd()
            return f(*a, **kw)
        finally:
            os.chdir(cwd)
    def decorator_update(dest, src):
        dest.__name__ = src.__name__
        dest.__doc__ = src.__doc__
        dest.__dict__.update(src.__dict__)
    decorator_update(cwd_preserver, f)
    return cwd_preserver

def selfupdate(onlyargs=False, exclude=[]):
    """Call in any method to set instance attributes from local variables.
    
    @param onlyargs: If True, use only named arguments
    
    @param exclude: Names of other variables to exclude.

    @note: Instance to update is assumed to be first argument of the caller.
    
    @note: Equivalent to vars(self).update(vars()); del self.self
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

def usetempfile(function):
    """Unittest method decorator.  Decorates a method of the
    form f(self, fname) so that fname is replaced with a temporary
    file which is deleted at the end of the call. 
    """
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


