#!/usr/bin/python
#!/export/home/medscan/local32/bin/python

"""MScanner XML-RPC service

Query medline with a corpus of citations.

"""

__docformat__ = "restructuredtext en"
from datetime import datetime
from DocXMLRPCServer import DocCGIXMLRPCRequestHandler
import os
from path import path
import re
import subprocess
import sys
import time
import unittest
from xmlrpclib import Fault

mtree = True

if mtree:
    source = path("/export/home/medscan/source")
    output = path("/export/apps/medline/htdocs/output")
    python = path("/export/home/medscan/local32/bin/python")
else:
    source = path("C:/Documents and Settings/Graham/My Documents/data/mscanner/source")
    output = path("C:/cygwin/srv/www/htdocs/mscanner/output")
    python = "C:/Python24/python.exe"
    
bin = source / "bin"
pidfile = source.parent / "data" / "cache" / "mscanner.pid"
lastlog = path("/tmp/mscanner-lastlog.txt")

def check_batchid(batchid):
    """
    Check that a given string is of the YYYYmmdd-HHMMSS form

    :param batchid: String for to check

    :raise Fault(2): Batch ID is invalid
    """
    if not re.match(r"^[A-Za-z0-9-]+$", batchid):
        raise Fault(2, "Invalid batch id %s" % batchid)

def readpid():
    """
    Get status from the pidfile of the MScanner spawned process

    :return: PID, progress, total, batch ID, start time
    """
    if not pidfile.isfile():
        return None, None, None, None, None
    lines = pidfile.lines()
    return int(lines[0]), int(lines[1]), int(lines[2]), lines[3].strip(), float(lines[4])

def generate_batch(positives):
    """
    Creates a batch directory and returns the batch ID

    :param positives: Text to be written to positives.txt

    :raise Fault(2): Directory for the batch already exists

    :return: The new batch id
    """
    # Create batch dir
    batchid = datetime.now().strftime("%Y%m%d-%H%M%S")
    bdir = output / batchid
    if bdir.exists():
        raise Fault(2, "Output for %s already exists" % batchid)
    bdir.mkdir()
    # Output positive PMIDs (catch errors in script)
    (bdir / "positives.txt").write_lines([str(p) for p in positives])
    return batchid

class MScannerService:

    def getStatus(self, batchid):
        """
        Returns a mapping with various optional status keys

        :param batchid: A batch ID or an empty string

        :return: A key:value mapping

        Result Keys
        ===========

        - `batchid`: Echo of the batchid parameter.

        - `lastlog`: String with log contents (given if status is 'notfound')

        - `status`: May be 'done', 'busy', 'free' or 'notfound'.
        'done' indicates the batch is complete.  'busy' indicates that
        it is in progress.  'free' is returned if an empty string is
        given and no batch is in progress.  'notfound' indicates if no
        batch with the given ID exists.  The fields below are only
        returned if the status is 'busy':

        - `elapsed`: Float for seconds elapsed since start of the batch.

        - `type`: Either 'query' or 'validation'. 

        - `progress`: Integer for number of steps completed

        - `total`: Integer for total number of steps required

        """
        if batchid != "":
            check_batchid(batchid)
        pid, progress, total, dummy_batchid, started = readpid()
        idx = output / batchid / "index.html"
        if idx.isfile():
            # Batch is done
            return dict(
                batchid=batchid,
                status="done"
                )
        else:
            if pid is not None:
                # Busy with this batch
                return dict(
                    batchid=dummy_batchid,
                    status="busy",
                    elapsed=time.time()-started,
                    progress=progress,
                    total=total,
                    )
            elif batchid == "":
                # No pidfile, no batch = free time
                return dict(
                    batchid=batchid,
                    status="free",
                    )
            else:
                # Didn't find the batch
                result = dict(
                    batchid=batchid,
                    status="notfound",
                    )
                if lastlog.isfile():
                    result["lastlog"] = lastlog.text()
                return result

    def listBatches(self):
        """Get IDs of completed batches

        :return: A list of batch IDs
        """
        result = []
        for d in output.dirs():
            try:
                check_batchid(d.basename())
                # Delete broken batches (nothing is running)
                if (d / "index.html").isfile():
                    result.append(str(d.basename()))
                else:
                    if not pidfile.isfile():
                        d.rmtree()
            except Fault:
                continue
        return result

    def deleteBatch(self, batchid):
        """Delete a batch

        :param batchid: ID of the batch to delete

        :raise Fault(3): Batch ID is invalid
        """
        check_batchid(batchid)
        todel = output / batchid
        if not todel.isdir():
            raise Fault(2, "Directory %s does not exist" % todel)
        todel.rmtree(todel)

    def query(self, positives, pseudocount, limit, threshold):
        """Query Medline with a corpus of citations

        :param positives: List of integer PubMed IDs

        :param pseudocount: Float with the Bayesian pseudocount vlaue
        to use (betwen 0.0 and 1.0)

        :param limit: Maximum number of results to return (between 1
        and 10000).

        :param threshold: Score threshold for results (may result in
        fewer results than specified in limit).

        :raise Fault(2): One or more parameters are invalid

        :raise Fault(3): The scanner service is busy

        :return: Batch ID for the query
        """
        if pidfile.isfile():
            raise Fault(3, "Scanner is busy")
        # Validate numeric arguments
        if not isinstance(positives, list):
            raise Fault(2, "Positives should be a list")
        for p in positives:
            if not isinstance(p, int):
                raise Fault(2, "Positives should be a list of integers")
        if not isinstance(limit, int) or limit < 1 or limit > 10000:
            raise Fault(2, "Limit is < 1 or > 10000")
        if not isinstance(pseudocount, float) or pseudocount < 0.0 or pseudocount > 1.0:
            raise Fault(2, "Pseudocount is < 0.0 or > 1.0")
        if not isinstance(threshold, float) and not isinstance(threshold, int):
            raise Fault(2, "Threshold needs to be a number")
        # Run query
        batchid = generate_batch(positives)
        try:
            lastlog.remove()
        except:
            pass
        subprocess.Popen(
            [python, bin / "query.py", batchid, str(pseudocount), str(limit), str(threshold)],
            stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
        lastlog.chmod(0666)
        return batchid

    def validate(self, positives, negatives, nfolds, pseudocount):
        """Validate query corpus for performance statistics and tuning

        :param positives: List of integer PubMed IDs

        :param negatives: Integer for the number of random PubMed IDs
        to use as a negative corpus (between 10 and 500000).
        
        :param nfolds: Number of validation folds to use (between 2
        and 10).

        :param pseudocount: Float with the Bayesian pseudocount vlaue
        to use (betwen 0.0 and 1.0)

        :raise Fault(2): One or more parameters are invalid

        :raise Fault(3): If the scanner service is busy

        :return: Batch ID for the validation task
        """
        if pidfile.isfile():
            raise Fault(3, "Scanner is busy")
        # Validate numeric arguments
        if not isinstance(positives, list):
            raise Fault(2, "Positives should be a list")
        for p in positives:
            if not isinstance(p, int):
                raise Fault(2, "Positives should be a list of integers")
        if not isinstance(negatives, int) or negatives < 10 or negatives > 500000:
            raise Fault(2, "Negatives count is < 10 or > 500000") 
        if not isinstance(nfolds, int) or nfolds < 2 or nfolds > 10:
            raise Fault(2, "Number of folds is < 2 or > 10")
        if not isinstance(pseudocount, float) or pseudocount < 0.0 or pseudocount > 1.0:
            raise Fault(2, "Pseudocount is < 0.0 or > 1.0")
        # Run validation
        batchid = generate_batch(positives)
        try:
            lastlog.remove()
        except:
            pass
        subprocess.Popen(
            [python, bin / "validate.py", batchid, str(negatives), str(nfolds), str(pseudocount)],
            stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
        lastlog.chmod(0666)
        return batchid

class MedscanTests(unittest.TestCase):
    def test(self):
        s = MScannerService()
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
    handler.set_server_title("MScanner")
    handler.set_server_name("MScanner: Query medline using a corpus of citations")
    handler.set_server_documentation(__doc__)
    handler.register_instance(MScannerService())
    handler.handle_request()
