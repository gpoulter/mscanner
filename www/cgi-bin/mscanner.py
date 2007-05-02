#!/usr/bin/python
#!/export/home/medscan/local32/bin/python

"""MScanner XML-RPC service

Query medline with a corpus of citations.

"""

__docformat__ = "restructuredtext en"
from datetime import datetime
from DocXMLRPCServer import DocCGIXMLRPCRequestHandler
import md5
import os
from path import path
import re
import subprocess
import sys
import time
import unittest
from xmlrpclib import Fault

mtree = False

if mtree:
    source = path("/export/home/medscan/source")
    mailer = source.parent / "data" / "cache" / "emails.txt"
    output = path("/export/apps/medline/htdocs/output")
    python = path("/export/home/medscan/local32/bin/python")
    lastlog = path("/tmp/mscanner-lastlog.txt")
else:
    source = path("C:/Documents and Settings/Graham/My Documents/data/mscanner/source")
    mailer = path("C:/Documents and Settings/Graham/My Documents/data/mscanner/data/cache")
    output = path("C:/Documents and Settings/Graham/My Documents/data/mscanner/data/output")
    lastlog = path("C:/cygwin/tmp/mscanner-lastlog.txt")
    python = "C:/Python24/python.exe"
    
bin = source / "bin"
pidfile = source.parent / "data" / "cache" / "mscanner.pid"

def check_batchid(batchid):
    """
    Check that a given string batch ID matches the batch ID regexp

    :param batchid: String to check

    :raise Fault(2): Batch ID is invalid
    """
    if not re.match(r"^[A-Za-z0-9_.-]+$", batchid):
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

def initialise_batch_dir(batchid, delcode, positives):
    """
    Creates a batch directory and returns the batch ID

    :param positives: Text to be written to positives.txt

    :raise Fault(2): Directory for the batch already exists

    :return: The new batch id
    """
    # Create batch dir
    #batchid = datetime.now().strftime("%Y%m%d-%H%M%S")
    bdir = output / batchid
    if bdir.exists():
        raise Fault(2, "A batch with title %s already exists" % batchid)
    bdir.mkdir()
    # Output positive PMIDs (catch errors in script)
    (bdir / "positives.txt").write_lines([str(p) for p in positives])
    # Write hash of deletion code (checked in deleteBatch)
    (bdir / "delcode").write_bytes(md5.new(delcode).hexdigest())

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
                # Delete broken batches which are not in-progress
                if (d / "index.html").isfile():
                    result.append(str(d.basename()))
                elif not pidfile.isfile():
                    d.rmtree()
            except Fault:
                continue
        return result

    def addMailer(self, email):
        """Queue an e-mail alert when scanner is available

        :param email: E-mail address to be alerted
        """
        mailer.write_lines(email, append=True)

    def deleteBatch(self, batchid, delcode):
        """Delete a batch

        :param batchid: ID of the batch to delete

        :param delcode: Deletion code of the batch (fails if incorrect)

        :raise Fault(3): Batch ID is invalid
        """
        check_batchid(batchid)
        dir_to_rm = output / batchid
        if not dir_to_rm.isdir():
            raise Fault(2, "Directory %s does not exist" % dir_to_rm)
        if md5.new(delcode).hexdigest() != (dir_to_rm / "delcode").bytes():
            raise Fault(2, "Deletion code %s is incorrect" % delcode)
        todel.rmtree(dir_to_rm)

    def query(self, batchid, delcode, positives, pseudocount, limit, threshold):
        """Query Medline with a corpus of citations

        :param batchid: Title for the batch (used as directory name)

        :param delcode: Deletion code for the batch (stop malicious deletions)

        :param positives: List of integer PubMed IDs

        :param pseudocount: Float with the Bayesian pseudocount vlaue
        to use (betwen 0.0 and 1.0)

        :param limit: Maximum number of results to return (between 1
        and 10000).

        :param threshold: Score threshold for results (may result in
        fewer results than specified in limit).

        :raise Fault(2): One or more parameters are invalid

        :raise Fault(3): The scanner service is busy
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
        initialise_batch_dir(batchid, delcode, positives)
        try:
            lastlog.remove()
        except:
            pass
        subprocess.Popen(
            [python, bin / "query.py", batchid, str(pseudocount), str(limit), str(threshold)],
            stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
        lastlog.chmod(0666)

    def validate(self, batchid, delcode, positives, negatives, nfolds, pseudocount, alpha):
        """Validate query corpus for performance statistics and tuning

        :param batchid: Title for the batch (used as directory name)

        :param delcode: Deletion code for the batch (stop malicious deletions)

        :param positives: List of integer PubMed IDs

        :param negatives: Integer for the number of random PubMed IDs
        to use as a negative corpus (between 10 and 500000).
        
        :param nfolds: Number of validation folds to use (between 2
        and 10).

        :param pseudocount: Bayesian pseudocount (betwen 0.0 and 1.0)

        :param alpha: Alpha-parameter for F-Measure tradeoff (between
        0.0 and 1.0, where alpha=0.5 maximises standard F-Measure)

        :raise Fault(2): One or more parameters are invalid

        :raise Fault(3): If the scanner service is busy
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
        if not isinstance(alpha, float) or alpha < 0.0 or alpha > 1.0:
            raise Fault(2, "F-Measure alpha is < 0.0 or > 1.0")
        # Run validation
        initialise_batch_dir(batchid, delcode, positives)
        try:
            lastlog.remove()
        except:
            pass
        subprocess.Popen(
            [python, bin / "validate.py", batchid, str(negatives), str(nfolds), str(pseudocount), str(alpha)],
            stdout=lastlog.open("w"), stderr=subprocess.STDOUT)
        lastlog.chmod(0666)

class MedscanTests(unittest.TestCase):
    def test(self):
        s = MScannerService()
        #self.assertEqual(listBatches(), [])
        self.assertRaises(Fault, s.deleteBatch, "noexist", "abc")
        self.assertEqual(s.getStatus("11223344-112233"), dict(batchid="11223344-112233", status="notfound"))
        # Test query
        s.query("test-query", "abc", [11809184,12069159,9744524,1960624], pseudocount=0.1, limit=5, threshold=0)
        time.sleep(2)
        self.assertEqual(s.getStatus("test-query")["status"], "busy")
        id, status = os.waitpid(int(pidfile.lines()[0].strip()),0)
        self.assertEqual(s.getStatus("test-query"), dict(batchid="test-query", status="done"))
        self.assertEqual(s.getStatus(""), dict(status="free"))
        self.assertRaises(Fault, s.deleteBatch, "test-query", "def")
        deleteBatch("test-query", "abc")
        # Test validation
        batchid = s.validate("test-valid", "def", [11809184,12069159,9744524,1960624],
                             negatives=50, nfolds=5, pseudocount=0.1, alpha=0.5)
        time.sleep(1)
        self.assertEqual(s.getStatus(batchid)["status"], "busy")
        id, status = os.waitpid(int(pidfile.lines()[0].strip()),0)
        self.assertEqual(s.getStatus(batchid), dict(batchid=batchid, status="done"))
        self.assertEqual(s.getStatus(""), dict(status="free"))
        deleteBatch("test-valid", "def")

if __name__ == "__main__":
    unittest.main()
    sys.exit(0)
    handler= DocCGIXMLRPCRequestHandler()
    handler.set_server_title("MScanner")
    handler.set_server_name("MScanner: Query medline using a corpus of citations")
    handler.set_server_documentation(__doc__)
    handler.register_instance(MScannerService())
    handler.handle_request()
