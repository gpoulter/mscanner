#!/usr/bin/env python

"""
Web.Py controller for the MScanner web interface (view is the
template code, and model is the command line scripts that are called).

                                   
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

import web
from web.utils import Storage
import re
import time
import platform

from templates import front, form, output, status

# In case of windows-specific handling
windows = (platform.system() == "Windows")

# Set error handler to be informative
web.webapi.internalerror = web.debugerror

# Configure URLs for the application
urls = (
    '/', 'web_front',
    '/form', 'web_form',
    '/status', 'web_status',
    '/output/(.*)', 'web_output',
)

defaults = dict(
    baseurl = "",
    inc = lambda x: "templates/" + x,
    )

form_defaults = Storage(
    dataset = "",
    code = "",
    pseudocount = 0.01,
    positives = "",
    threshold = 0,
    limit = 500,
    nfold = 10,
    negatives = 50000,
    alpha = 0.5
    )

class web_front:
    """Front page of the site"""
    def GET(self):
        """Return the front page, with a status box if there is something
        currently running"""
        t = front.front(searchList=defaults)
        print t

class web_form:
    """Main submission form for queries or validation"""
    def GET(self):
        """Return the form with default values"""
        t = form.form(searchList=[defaults, form_defaults])
        t.status = None
        t.errors = None
        print t
        
    def POST(self):
        """Submit the form.  If errors, return with the same values
        but errors listed at the top.
        """
        # Get form values for fields listed in form_defaults
        i = web.input(*form_defaults.keys())
        errors = self.validateForm(i)
        t = form.form(searchList=[defaults, testdata.form])
        t.status = testdata.status
        print t
        
    def validateForm(self, input):
        """Given a Storage object with keys as in form_defaults,
        return None if there are no errors, or a list of strings
        describing the errors"""
        errors = []
        if re.match(r"^[a-z._-]+$", input.dataset) is None:
            errors.append("Data set code does not match the [a-z._-]+"+ 
            "restriction (due to it being used as the directory name).")
       

class web_status:
    """Status page"""
    def GET(self):
        t = status.status(searchList=defaults)
        t.status = testdata.status
        t.output = "BLAH BLAH BLAH"
        print t
        
    def POST(self):
        i = web.input("dataset", "code")
        t = status.status(searchList=defaults)
        t.status = testdata.status
        t.output = "BLAH BLAH BLAH"
        print t

class web_output:
    """Page with index of outputs"""
    def GET(self, name):
        t = output.output(searchList=defaults)
        t.datasets = testdata.datasets
        print t
        
    def POST(self):
        t = output.output(searchList=defaults)
        t.datasets = testdata.datasets
        print t

if __name__ == "__main__":
    web.run(urls, globals())
