#!/usr/bin/env python

"""CGI XML-RPC medscanner service

Note that the script expects to be executed from its directory.

"""

import os
import re
import sys
import subprocess
import unittest
import time
from path import path
from datetime import datetime
from xmlrpclib import Fault
from DocXMLRPCServer import DocCGIXMLRPCRequestHandler

#mtree = True
mtree = False

if mtree:
    source = path("/export/home/medscan/source")
    www = path("/export/apps/medscan")
else:
    source = path("/c/Documents and Settings/Graham/My Documents/mutable/medscanner/source/")
    www = source / "www"
    
cgi_bin = www / "cgi-bin"
htdocs = www / "htdocs"
output_dir = htdocs / "output"
bin_dir = source / "bin"
pidfile = path("/var/run/medscan.pid")
lastlog = cgi_bin / "lastlog.txt"

def check_batchid(batchid):
    """Batch IDs must be of YYYYmmdd-HHMMSS form"""
    if not re.match(r"^\d{8}-\d{6}$", batchid):
        raise Fault(2, "Invalid batch id %s" % batchid)

def readpid():
    """Return pid, batchid, start_time from pidfile"""
    if not pidfile.isfile():
        return None, None, None
    lines = pidfile.lines()
    return int(lines[0]), lines[1].strip(), float(lines[2])

def getStatus(batchid):
    """Returns a dict which has at least the status element, which may
    be 'done', 'busy', 'free', or 'notfound'.  All except 'free' (only
    returned when batchid is empty string) send back a batchid
    element.  'busy' includes a started element in seconds since
    epoch.
    """
    if batchid != "": check_batchid(batchid)
    pid, pbatch, started = readpid()
    idx = output_dir / batchid / "index.html"
    if idx.isfile():
        return dict(batchid=batchid, status="done")
    else:
        if pid is not None:
            return dict(batchid=pbatch, status="busy", started=started)
        elif batchid == "":
            return dict(status="free")
        else:
            result = dict(batchid=batchid, status="notfound")
            if lastlog.isfile(): result["lastlog"] = lastlog.text()
            return result

def listBatches():
    """Return list of batch IDs whose results are available"""
    if not pidfile.isfile():
        for d in output_dir.dirs():
            check_batchid(d.basename())
            if not (d / "index.html").isfile():
                d.rmtree()
    return [ str(d.basename()) for d in output_dir.dirs() ]

def deleteBatch(batchid):
    """Delete output corresponding to given batchid.  Fault 2 if
    directory does not exist.
    """
    check_batchid(batchid)
    todel = output_dir / batchid
    if not todel.isdir():
        raise Fault(2, "Directory %s does not exist" % todel)
    todel.rmtree(todel)

def query(positives, pseudocount, limit):
    """Begin query process, returning batchid"""
    # Validate numeric arguments
    if not isinstance(limit, int) or limit < 1 or limit > 10000:
        raise Fault(2, "Limit is < 1 or > 10000")
    if not isinstance(pseudocount, float) or pseudocount < 0.0 or pseudocount > 1.0:
        raise Fault(2, "Pseudocount is < 0.0 or > 1.0")
    # Create batch dir
    batchid = datetime.now().strftime("%Y%m%d-%H%M%S")
    bdir = output_dir / batchid
    if bdir.exists(): raise Fault(2, "Output for %s already exists" % batchid)
    bdir.mkdir()
    # Output positive PMIDs (catch errors in script)
    (bdir / "positives.txt").write_text(positives)
    batchid = datetime.now().strftime("%Y%m%d-%H%M%S")
    # Run query
    subprocess.Popen(
        ["python", bin_dir / "query.py", batchid, str(pseudocount), str(limit)],
        stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
    return batchid

def validate(positives, negatives, nfolds, recall, pseudocount):
    """Begin validation process, returning batchid"""
    # Validate numeric arguments
    if not isinstance(negatives,int) or negatives < 10 or negatives > 500000:
        raise Fault(2, "Negatives count is < 10 or > 500000") 
    if not isinstance(nfolds, int) or nfolds < 2 or nfolds > 10:
        raise Fault(2, "Number of folds is < 2 or > 10")
    if not isinstance(recall, float) or recall < 0.0 or recall > 1.0:
        raise Fault(2, "Recall is < 0.0 or > 1.0")
    if not isinstance(pseudocount, float) or pseudocount < 0.0 or pseudocount > 1.0:
        raise Fault(2, "Pseudocount is < 0.0 or > 1.0")
    # Create batch dir
    batchid = datetime.now().strftime("%Y%m%d-%H%M%S")
    bdir = output_dir / batchid
    if bdir.exists(): raise Fault(2, "Output for %s already exists" % batchid)
    bdir.mkdir()
    # Output positive PMIDs
    (bdir / "positives.txt").write_text(positives)
    # Run query
    subprocess.Popen(
        ["python", bin_dir / "validate.py", batchid, str(negatives), str(nfolds), str(recall), str(pseudocount)],
        stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
    return batchid

class MedscanTests(unittest.TestCase):
    def test(self):
        #self.assertEqual(listBatches(), [])
        self.assertRaises(Fault, deleteBatch, "noexist")
        self.assertRaises(Fault, deleteBatch, "11223344-112233")
        self.assertEqual(getStatus("11223344-112233"), dict(batchid="11223344-112233", status="notfound"))
        # Test query
        batchid = query("11809184\n12069159\n9744524\n1960624\n", 0.1, 5)
        time.sleep(2)
        self.assertEqual(getStatus(batchid)["status"], "busy")
        id, status = os.waitpid(int(pidfile.lines()[0].strip()),0)
        self.assertEqual(getStatus(batchid), dict(batchid=batchid, status="done"))
        self.assertEqual(getStatus(""), dict(status="free"))
        deleteBatch(batchid)
        # Test validation
        batchid = validate("11809184\n12069159\n9744524\n1960624\n", 50, 5, 0.9, 0.1)
        time.sleep(1)
        self.assertEqual(getStatus(batchid)["status"], "busy")
        id, status = os.waitpid(int(pidfile.lines()[0].strip()),0)
        self.assertEqual(getStatus(batchid), dict(batchid=batchid, status="done"))
        self.assertEqual(getStatus(""), dict(status="free"))
        deleteBatch(batchid)

if __name__ == "__main__":
    #unittest.main()
    handler= DocCGIXMLRPCRequestHandler()
    handler.set_server_title("WebTwain")
    handler.set_server_name("WebTwain - scan documents over the web")
    handler.set_server_documentation("This is an automatically generated description of the scanner's XML-RPC interface")
    handler.register_function(getStatus)
    handler.register_function(listBatches)
    handler.register_function(deleteBatch)
    handler.register_function(validate)
    handler.register_function(query)
    handler.handle_request()
