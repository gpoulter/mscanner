"""Graham's Cheetah frontend"""

                                     
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

import warnings
warnings.simplefilter("ignore", UserWarning)

import logging as log
from path import path
import sys
import Cheetah.Template

class FileTransaction(file):
    """Transaction for Cheetah templates to output direct-to-file.
    
    Cheetah defaults to DummyTransaction which creates a huge list and
    joins them up to create a string.  This is way slower than writing to
    file directly.
    
    Usage::
        ft = FileTransaction("something.html","wb")
        Template().respond(ft)
        ft.close()
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



class TemplateMapper:
    """Attribute-style access to non-compiled Cheetah templates..
    
    This hides the distinction between importing a pre-compiled template from a
    module versus parsing a raw file template in a subdirectory. Seemed like a
    good idea at the time."""

    def __init__(self, 
                 root=path("templates"), 
                 ext=".tmpl", 
                 module=None, 
                 gvars={}, 
                 template_class=Cheetah.Template.Template, 
                 kwargs={}):
        """Initialise template mapper.
        
        @param root: Path to root template directory
        
        @param ext: Template extension
        
        @param module: Root modulue for templates (root is then not used)
        
        @param gvars: Global variables dictionary
        
        @param template_class: Defaults to Cheetah.Template.Template
        
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
        """Attribute access to mapper creates a callable template object"""
        # Instantiate using self.root / name
        if self.module is None:
            return self.CallableTemplate(
                template = self.root / name + self.ext, 
                default_search = self.gvars,
                template_class = self.template_class)
        # Instantiate using self.module.name.name
        else:
            return self.CallableTemplate(
                template = getattr(getattr(self.module, name), name), 
                default_search = self.gvars,
                template_class = self.template_class)


    class CallableTemplate:
        """Proxy for a Cheetah template.
        
        Calling this proxy template fills it, with correct constructor
        invocation for file, string or compiled templates."""
        
        def __init__(self, 
                     template, 
                     default_search={}, 
                     template_class=Cheetah.Template.Template, 
                     kwargs={}):
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
            
            @param pvars: Dictionaries to add to the searchlist
            @param kvars: Keywords to add to search list."""
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
