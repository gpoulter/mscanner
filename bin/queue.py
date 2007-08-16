#!/usr/bin/env python

"""Queueing facility for the web frontend

Web frontend creates a YYYYMMDDHHMMSS.txt descriptor file in the queue
directory, formatted as follows for query:

#operation = query
#dataset = Whatever
#limit = 500
#threshold = 10.3
#timestamp = 23424123.3
804133
3214241
...

And for cross validation:

#operation = validate
#dataset = Whatever
#nfolds = 10
#numnegs = 100000
#alpha = 0.5
#timestamp = 23424123.3
804133
3214241
...

Queueing program checks this directory every second for new files, and starts a
query or validation operation. When the operation completes, the descriptor
file is moved to the output.
 
                                   
"""

__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""

import logging as log
import os
from path import path
import sys
import time

from mscanner.configuration import rc, initLogger
from mscanner.validenv import ValidationEnvironment
from mscanner.queryenv import QueryEnvironment
from mscanner.scorefile import readDescriptor, writeDescriptor

rc.query_report_dir = lambda: rc.web_report_dir / rc.dataset
rc.valid_report_dir = lambda: rc.web_report_dir / rc.dataset

def mainloop():
    """Look for descriptor files, and perform query or validation
    as requested"""
    qenv = QueryEnvironment()
    venv = ValidationEnvironment()
    # Initialise the article list for later
    venv.article_list
    while True:
        listing = rc.queue_path.files("*.txt")
        if listing:
            descriptor = listing[0]
            rc.update(readDescriptor(descriptor))
            log.info("Performing %s on descriptor %s", 
                     rc.operation, descriptor.basename())
            try:
                # REMOVE SLEEP FOR MAPLES - TESTING ONLY
                time.sleep(20)
                if rc.operation == "query":
                    qenv.standardQuery(descriptor)
                    descriptor.move(rc.query_report_dir / "descriptor.txt")
                elif rc.operation == "validate":
                    venv.standardValidation(descriptor)
                    descriptor.move(rc.valid_report_dir / "descriptor.txt")
            except ValueError, e:
                log.error(e)
                descriptor.remove()
        time.sleep(1)

def makeTestQueue():
    from mscanner.scorefile import readPMIDs
    from mscanner.support.storage import Storage
    pmids = list(readPMIDs(rc.corpora / "genedrug-small.txt"))
    ts = lambda x, d: rc.queue_path / \
    (time.strftime("%Y%m%d%H%M%S", time.gmtime(x)) + "-" + d + ".txt")
    params = Storage(
        operation = "validate", dataset = "gdsmallvalidq",
        nfolds = 5, numnegs = 1000, alpha = 0.6,
        timestamp = time.time(), limit = 500, threshold = 0)
    writeDescriptor(ts(params.timestamp, params.dataset), pmids, params)
    params.timestamp += 5
    params.operation = "query"
    params.dataset = "gdsmallqueryq"
    writeDescriptor(ts(params.timestamp, params.dataset), pmids, params)

if __name__ == "__main__":
    initLogger()
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        makeTestQueue()
    mainloop()