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
    output = path("/export/apps/medscan/htdocs/output")
else:
    source = path("/home/graham/data/mscanner/source")
    output = path("/srv/www/htdocs/mscanner/output")
    
bin = source / "bin"
pidfile = path("/var/run/mscanner.pid")
lastlog = path("/tmp/mscanner-lastlog.txt")

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

def generate_batch(positives):
    # Create batch dir
    batchid = datetime.now().strftime("%Y%m%d-%H%M%S")
    bdir = output / batchid
    if bdir.exists(): raise Fault(2, "Output for %s already exists" % batchid)
    bdir.mkdir()
    # Output positive PMIDs (catch errors in script)
    (bdir / "positives.txt").write_text(positives)
    return batchid

class MScannerService:

    def getStatus(self, batchid):
        """Returns a dict which has at least the status element, which may
        be 'done', 'busy', 'free', or 'notfound'.  All except 'free' (only
        returned when batchid is empty string) send back a batchid
        element.  'busy' includes a started element in seconds since
        epoch.
        """
        if batchid != "": check_batchid(batchid)
        pid, pbatch, started = readpid()
        idx = output / batchid / "index.html"
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

    def listBatches(self):
        """Return list of batch IDs whose results are available"""
        result = []
        if not pidfile.isfile():
            for d in output.dirs():
                try:
                    check_batchid(d.basename())
                except Fault:
                    continue
                if not (d / "index.html").isfile():
                    d.rmtree()
                else:
                    result.append(str(d.basename()))
        return result

    def deleteBatch(self, batchid):
        """Delete output corresponding to given batchid.  Fault 2 if
        directory does not exist.
        """
        check_batchid(batchid)
        todel = output / batchid
        if not todel.isdir():
            raise Fault(2, "Directory %s does not exist" % todel)
        todel.rmtree(todel)

    def query(self, positives, pseudocount, limit, threshold):
        """Begin query process, returning batchid"""
        # Validate numeric arguments
        if not isinstance(limit, int) or limit < 1 or limit > 10000:
            raise Fault(2, "Limit is < 1 or > 10000")
        if not isinstance(pseudocount, float) or pseudocount < 0.0 or pseudocount > 1.0:
            raise Fault(2, "Pseudocount is < 0.0 or > 1.0")
        if not isinstance(threshold, float) and not isinstance(threshold, int):
            raise Fault(2, "Threshold needs to be a number")
        # Run query
        batchid = generate_batch(positives)
        subprocess.Popen(
            ["python", bin / "query.py", batchid, str(pseudocount), str(limit), str(threshold)],
            stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
        return batchid

    def validate(self, positives, negatives, nfolds, pseudocount):
        """Begin validation process, returning batchid"""
        # Validate numeric arguments
        if not isinstance(negatives,int) or negatives < 10 or negatives > 500000:
            raise Fault(2, "Negatives count is < 10 or > 500000") 
        if not isinstance(nfolds, int) or nfolds < 2 or nfolds > 10:
            raise Fault(2, "Number of folds is < 2 or > 10")
        if not isinstance(pseudocount, float) or pseudocount < 0.0 or pseudocount > 1.0:
            raise Fault(2, "Pseudocount is < 0.0 or > 1.0")
        # Run validation
        batchid = generate_batch(positives)
        subprocess.Popen(
            ["python", bin / "validate.py", batchid, str(negatives), str(nfolds), str(pseudocount)],
            stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
        return batchid

class MedscanTests(unittest.TestCase):
    def test(self):
        s = MScannerSice()
        #self.assertEqual(listBatches(), [])
        self.assertRaises(Fault, s.deleteBatch, "noexist")
        self.assertRaises(Fault, s.deleteBatch, "11223344-112233")
        self.assertEqual(s.getStatus("11223344-112233"), dict(batchid="11223344-112233", status="notfound"))
        # Test query
        batchid = s.query("11809184\n12069159\n9744524\n1960624\n", 0.1, 5)
        time.sleep(2)
        self.assertEqual(s.getStatus(batchid)["status"], "busy")
        id, status = os.waitpid(int(pidfile.lines()[0].strip()),0)
        self.assertEqual(s.getStatus(batchid), dict(batchid=batchid, status="done"))
        self.assertEqual(s.getStatus(""), dict(status="free"))
        deleteBatch(batchid)
        # Test validation
        batchid = s.validate("11809184\n12069159\n9744524\n1960624\n", 50, 5, 0.9, 0.1)
        time.sleep(1)
        self.assertEqual(s.getStatus(batchid)["status"], "busy")
        id, status = os.waitpid(int(pidfile.lines()[0].strip()),0)
        self.assertEqual(s.getStatus(batchid), dict(batchid=batchid, status="done"))
        self.assertEqual(s.getStatus(""), dict(status="free"))
        deleteBatch(batchid)

if __name__ == "__main__":
    #unittest.main()
    handler= DocCGIXMLRPCRequestHandler()
    handler.set_server_title("WebTwain")
    handler.set_server_name("WebTwain - scan documents over the web")
    handler.set_server_documentation("This is an automatically generated description of the scanner's XML-RPC interface")
    handler.register_instance(MScannerService())
    handler.handle_request()
