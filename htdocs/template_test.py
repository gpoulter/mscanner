#!/usr/bin/env python

"""Web.Py controller for testing the templates of the MScanner web interface.
It implements dummy responders and has no dependencies on the main MScanner
code.

                                   
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
import time
import platform

web.webapi.internalerror = web.debugerror

from templates import front, form, output, status

urls = (
    '/', 'web_front',
    '/form', 'web_form',
    '/status', 'web_status',
    '/output/(.*)', 'web_output',
)

testdata = Storage(

    form = Storage(
        dataset = "",
        code = "",
        pseudocount = 0.01,
        positives = "",
        threshold = 0,
        limit = 1000,
        nfold = 10,
        negatives = 100000,
        alpha = 0.5
        ),

    status = Storage(
        dataset = "TESTING",
        timestamp = time.time()-300.0,
        progress = 10,
        total = 20
    ),

    errors = [ "Error A", "Error B" ],
    
    datasets = [ Storage(dataset="ABC", timestamp=time.time()),
                Storage(dataset="DEF", timestamp=time.time()-200) ]
)

defaults = dict(
    baseurl = "",
    inc = lambda x: "templates/" + x,
    )

class web_front:
    """Front page of the site"""
    def GET(self):
        t = front.front(searchList=defaults)
        t.status = testdata.status
        print t

class web_form:
    """Main submission form for queries or validation"""
    def GET(self):
        t = form.form(searchList=[defaults, testdata.form])
        t.status = None
        t.errors = None
        print t
    def POST(self):
        i = web.input("dataset", "code", "pseudocount", "positives",
                      "threshold", "limit", "nfold", "negatives", "alpha")
        t = form.form(searchList=[defaults, testdata.form])
        t.status = testdata.status
        print t

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
    #output().GET()
    web.run(urls, globals())
