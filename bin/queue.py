#!/usr/bin/env python

"""Queueing facility for the web frontend

Web frontend creates a YYYYMMDDHHMMSS.txt descriptor file in the queue
directory, formatted as follows for query:

#operation = query
#dataset = Whatever
#limit = 500
#threshold = 10.3
#timestamp = 23424123
804133
3214241
...

And for cross validation:

#operation = validate
#dataset = Whatever
#nfolds = 10
#numnegs = 100000
#alpha = 0.5
#timestamp = 23424123
804133
3214241
...

Queueing program checks this directory every second for new files, and starts a
query or validation operation. When the operation completes, the descriptor
file is moved to the output.
 
                                   
"""

__license__ = """ This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import logging as log
import os
from path import path
import sys
import time

from mscanner.configuration import rc, initLogger
from mscanner.validenv import ValidationEnvironment
from mscanner.queryenv import QueryEnvironment

def sendmail(smtpserver, email):
    """Read e-mail addresses from mail file and send them a message saying
    MScanner is available.

    @param mailer: Path object to file of e-mail addresses
    """
    import smtplib
    import logging
    server = smtplib.SMTP(smtpserver)
    server.set_debuglevel(0)
    fromaddr = "nobody@mscanner.stanford.edu"
    logging.debug("Sending availability alert to %s", email)
    msg = "From: %s\r\nTo: %s\r\n\r\n" % (fromaddr, email)
    msg += """
This is a once-off notification - no record of this email address has
been kept.
"""
    try:
        server.sendmail(fromaddr, email, msg)
    except Exception:
        logging.debug("Failed to send to %s", email)
    server.quit()
    mailer.remove()

def readDescriptor(fpath):
    """Reads the descriptor file, by setting rc parameters.

    @note: The same descriptor file is the input PubMed ID list,
    since the PubMed-ID reader ignores #-lines.
    """
    converter = dict(operation=str, dataset=str, limit=int, 
                     threshold=float, nfolds=int, 
                     numnegs=int, alpha=float, timestamp=float) 
    f = file(fpath, "r")
    line = f.readline()
    while line.startswith("#"):
        key, value = line[1:].split(" = ",1)
        value = converter[key](value.strip())
        rc[key] = value
        line = f.readline()
    f.close()
    
def writeDescriptor(fpath, pmids, params):
    """Write parameters and PubMed IDs to the descriptor file"""
    f = file(fpath, "w")
    for key, value in params.iteritems():
        f.write("#" + key + " = " + str(value) + "\n")
    fpath.write_lines([str(pmid) for pmid in pmids])
    f.close()

def mainloop():
    """Look for descriptor files, and perform query or validation
    as requested"""
    qenv = QueryEnvironment()
    venv = ValidationEnvironment()
    while True:
        listing = rc.queue_path.files("*.txt")
        if listing:
            descriptor = listing[0]
            readDescriptor(descriptor)
            log.info("Perfrming %s on descriptor %s", 
                     rc.operation, descriptor.basename())
            if rc.operation == "query":
                qenv.standardQuery(descriptor)
                descriptor.move(rc.query_report_dir)
            elif rc.operation == "validate":
                venv.standardValidation(descriptor)
                descriptor.move(rc.valid_report_dir)
        time.sleep(1)

def makeTestQueue():
    from mscanner.scorefile import readPMIDs
    timestamp = time.time()
    destfile = rc.queue_path / time.strftime("%Y%m%d-%H%M%S.txt",  
                                             time.gmtime(timestamp))
    pmids = list(readPMIDs(rc.corpora / "genedrug-small.txt"))
    writeDescriptor(destfile, pmids, dict(
        operation = "query", limit = 500, threshold = 0, 
        timestamp = timestamp, dataset = "gdsmallqueue"))
    timestamp += 5
    destfile = rc.queue_path / time.strftime("%Y%m%d-%H%M%S.txt",  
                                             time.gmtime(timestamp))
    writeDescriptor(destfile, pmids, dict(
        operation = "validate", dataset = "gdsmallvalidqueue",
        nfolds = 5, numnegs = 1000, timestamp = timestamp, alpha = 0.6))

if __name__ == "__main__":
    initLogger()
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        makeTestQueue()
    mainloop()