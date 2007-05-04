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
from Cheetah.Template import Template

class FileTransaction(file):
    """Transaction for template which outputs to file.  Cheetah
    default uses DummyTransaction which creates a huge list and joins
    them up, which is way slow for large files.

    Call the template instance's 'respond' method with a
    FileTransaction instead of converting the template instance to a
    string.
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

    def __init__(self, root=path("templates"), ext=".tmpl", module=None, gvars={}, template_cls=Template, kwargs={}):
        """Initialise template mapper.

        @param gvars: Global variables dictionary
        @param root: Path to template directory
        @param ext: Template extension
        @param module: Templates are modules under this module instead of re-compiling
        @param kwargs: Keyword arguments to pass to __init__ of the template class 
        """
        self.gvars = gvars
        self.gvars["_inc"] = lambda x: str(self.root / x)
        self.root = path(root)
        self.ext = ext
        self.module = module
        self.template_cls = template_cls
        self.kwargs = kwargs

    def __getattr__(self, name):
        """Use attribute access to obtain a callable template object"""
        if self.module is None:
            return self.CallableTemplate(self.root / name + self.ext, self.gvars, self.template_cls)
        else:
            return self.CallableTemplate(getattr(getattr(self.module, name), name), self.gvars, self.template_cls)

    class CallableTemplate:
        
        def __init__(self, template, gvars={}, template_cls=Template, kwargs={}):
            """Initialise callable template with set of global variables"""
            self.template = template
            self.gvars = gvars
            self.template_cls = template_cls
            self.kwargs = kwargs

        def __call__(self, *pvars, **kvars):
            """Call template with self.gvars for globals, pvars adding to
            the searchlist, and kvars forming a dictionary of parameters."""
            slist = list(pvars) + [kvars, self.gvars]
            if isinstance(self.template, path):
                return self.template_cls(file=str(self.template), searchList=slist, **self.kwargs)
            elif isinstance(self.template, basestring):
                return self.template_cls(self.template, searchList=slist, **self.kwargs)
            elif issubclass(self.template, Template):
                return self.template(searchList=slist, **self.kwargs)
            else:
                raise ValueError("Unrecognised template")
