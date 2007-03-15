#!/usr/bin/env python

"""
Web.Py controller for the MScanner web interface.  The templates are
the view, and the command-line scripts are the model.

                                   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
"""

import web
from web.utils import Storage
import time
import platform

from mscanner.utils import TemplateMapper

windows = platform.system() == "Windows"
web.webapi.internalerror = web.debugerror

urls = (
    '/', 'front',
    '/form', 'form',
    '/status', 'status',
    '/output/(.*)', 'output',
    '/hello/(.*)', 'hello'
)

testing = Storage(
    form = Storage(
    batch = "",
    code = "",
    pseudocount = 0.01,
    positives = "",
    threshold = 0,
    limit = 1000,
    nfold = 10,
    negatives = 100000,
    alpha = 0.5),

    status = Storage(
    batch = "TESTING",
    start = time.time()-300.0,
    progress = 10,
    total = 20
    ),

    errors = [ "Error A", "Error B" ],
    
    batches = [ Storage(batch="ABC", start=time.time()),
                Storage(batch="DEF", start=time.time()-200) ]
    )

r = TemplateMapper(dict(root=""))

class front:
    def GET(self):
        print r.page(r.front(status=testing.status))

class form:
    def GET(self):
        print r.page(r.form(testing.form, status=None, errors=None))
    def POST(self):
        i = web.input("batch", "code", "pseudocount", "positives",
                      "threshold", "limit", "nfold", "negatives", "alpha")
        
        print r.page(r.form(testing.form, status=testing.status))

class status:
    def GET(self):
        print r.page(r.status(status=testing.status, output="BLAH BLAH\nBLAH"))
    def POST(self):
        i = web.input("batch", "code")
        print r.page(r.status(status=testing.status, output="BLAH BLAH\nBLAH"))

class output:
    def GET(self, name):
        print r.page(r.output(batches=testing.batches))
        print name
    def POST(self):
        print r.page(r.output(batches=testing.batches))

if __name__ == "__main__":
    #output().GET()
    web.run(urls, globals())
