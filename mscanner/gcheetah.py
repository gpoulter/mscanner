"""Graham's Cheetah frontend

TemplateMapper -- Convenient frontend to Cheetah templating engine
FileTransaction -- File subclass to output cheetah stuff straight to disk

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import warnings
warnings.simplefilter("ignore", UserWarning)

import logging as log
from path import path
import sys
from Cheetah.Template import Template

class FileTransaction(file):
    """Transaction for template which outputs to file.  
    
    @note: Cheetah defaults to DummyTransaction which creates a 
    huge list and joins them up to output a string.  This is way slow for 
    large files (list append is O(N^2))

    @note: With an instance of Template, call respond() with a FileTransaction
    instance, so the template uses FileTransaction.writeln() directly. 
    """
    
    def __init__(self, *pargs, **kwargs):
        file.__init__(self, *pargs, **kwargs)
        self.response = self

    def writeln(self):
        self.write(txt)
        self.write('\n')

    def getvalue(self):
        return None

    def __call__(self):
        return self

class TemplateMapper:
    """Provides attribute-style access to Cheetah templates, for
    either importing template modules from a package or compiling
    files in the template directory."""

    def __init__(self, root=path("templates"), ext=".tmpl", module=None, gvars={}, template_class=Template, kwargs={}):
        """Initialise template mapper.
        
        @param gvars: Global variables dictionary
        
        @param root: Path to root template directory
        
        @param ext: Template extension
        
        @param module: Root modulue for templates (root is then not used)
        
        @param kwargs: Keyword args for tmpl_class.__init__()
        """
        self.root = path(root)
        self.gvars = gvars
        self.gvars["_inc"] = lambda x: str(self.root / x)
        self.gvars["_log"] = lambda x: sys.stderr.write(x)
        self.ext = ext
        self.module = module
        self.template_class = template_class
        self.kwargs = kwargs

    def __getattr__(self, name):
        """Attribute access creates a callable template object"""
        # Instantiate using self.root / name
        if self.module is None:
            return self.CallableTemplate(
                template=self.root / name + self.ext, 
                default_search=self.gvars,
                template_class=self.template_class)
        # Instantiate using self.module.name.name
        else:
            return self.CallableTemplate(
                template=getattr(getattr(self.module, name), name), 
                default_search=self.gvars,
                template_class=self.template_class)

    class CallableTemplate:
        
        def __init__(self, template, default_search={}, template_class=Template, kwargs={}):
            """Initialise callable template with set of global variables

            @param template: path, string, or Template subclass
            
            @param default_search: Default search dictionary
            
            @param template_class: What to instantiate if template is path/string
            
            @param kwargs: Keyword arguments for Template constructor
            """
            self.template = template
            self.default_search = default_search
            self.template_class = template_class
            self.kwargs = kwargs

        def __call__(self, *pvars, **kvars):
            """Call template with self.gvars for globals.
            
            @param *pvars: Dictionaries to add to the searchlist
            @param **kvars: Keywords to add to search list."""
            slist = list(pvars) + [kvars, self.default_search]
            if isinstance(self.template, path):
                return self.template_class(
                    file=str(self.template), searchList=slist, **self.kwargs)
            elif isinstance(self.template, basestring):
                return self.template_class(
                    source=self.template, searchList=slist, **self.kwargs)
            elif issubclass(self.template, Template):
                return self.template(
                    searchList=slist, **self.kwargs)
            else:
                raise ValueError("Unrecognised template")
