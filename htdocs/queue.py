#!/usr/bin/env python

"""Queueing facility for the web frontend

Queueing program checks queue directory every second for new descriptor files,
and starts a query or validation operation. When the operation completes, the
descriptor file is moved to the output.

Basic file format for query::

    #operation = query
    #dataset = Whatever
    #limit = 500
    #threshold = 10.3
    #timestamp = 23424123.3
    804133
    3214241
    ...
    
Basic file format for cross validation::

    #operation = validate
    #dataset = Whatever
    #nfolds = 10
    #numnegs = 100000
    #alpha = 0.5
    #timestamp = 23424123.3
    804133
    3214241
    ...
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

import logging as log
import os
from path import path
import sys
import time

from mscanner.configuration import rc, initLogger
from mscanner import scorefile
from mscanner import validenv
from mscanner import queryenv

rc.query_report_dir = lambda: rc.web_report_dir / rc.dataset
rc.valid_report_dir = lambda: rc.web_report_dir / rc.dataset


def delete_old_outputs():
    """Find the oldest outputs in the web report directory and remove them"""
    max_keep = 100 # Max number of outputs to keep
    # Read outupt descriptors and sort by date
    output = [] # List of (timestamp, path)
    for fpath in rc.web_report_dir.dirs():
        if (fpath / rc.report_index).exists():
            descriptor = scorefile.readDescriptor(fpath / rc.report_descriptor)
            output.append((descriptor.timestamp, fpath))
    output.sort(reverse=True)
    # If more than max_keep in the list, delete the oldest ones
    for timestamp, outpath in output[max_keep:]:
        try:
            log.info("Cleaning out result %s", outpath.basename())
            for fname in outpath.files():
                fname.remove()
            outpath.rmdir()
        except OSError:
            # Whoops - failed to delete an output
            pass


def mainloop():
    """Look for descriptor files every second"""
    qenv = queryenv.QueryEnvironment()
    venv = validenv.ValidationEnvironment()
    venv.article_list # Initialise the article list for later
    loop_count = 0
    while True:
        listing = rc.queue_path.files("*.txt")
        if listing:
            descriptor = listing[0]
            rc.update(scorefile.readDescriptor(descriptor))
            log.info("Performing %s on descriptor %s", 
                     rc.operation, descriptor.basename())
            try:
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
        loop_count += 1
        if loop_count % 3600 == 0:
            # Every hour go delete old outputs
            delete_old_outputs()


def makeTestQueue():
    from mscanner.support.storage import Storage
    pmids = list(scorefile.readPMIDs(rc.corpora / "genedrug-small.txt"))
    ts = lambda x, d: rc.queue_path / \
    (time.strftime("%Y%m%d%H%M%S", time.gmtime(x)) + "-" + d + ".txt")
    params = Storage(
        operation = "validate", dataset = "gdsmallvalidq",
        nfolds = 5, numnegs = 1000, alpha = 0.6,
        timestamp = time.time(), limit = 500, threshold = 0)
    scorefile.writeDescriptor(ts(params.timestamp, params.dataset), pmids, params)
    params.timestamp += 5
    params.operation = "query"
    params.dataset = "gdsmallqueryq"
    scorefile.writeDescriptor(ts(params.timestamp, params.dataset), pmids, params)


if __name__ == "__main__":
    initLogger()
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        makeTestQueue()
    mainloop()