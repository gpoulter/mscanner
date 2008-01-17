#!/usr/bin/env python

"""Utility functions for MScanner database files and files containing lists of
PubMed IDs.

Usage::
    ./dbhelper.py function_name arg1 arg2 [...]

Please see the source code for the different functions.

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from __future__ import with_statement
from bsddb import db
from contextlib import closing
import numpy as nx
from path import path
import random
import sys

from mscanner.core import iofuncs
from mscanner.medline.FeatureDatabase import FeatureDatabase
from mscanner.medline.FeatureStream import FeatureStream, DateAsInteger
from mscanner.medline.Updater import Updater
from mscanner.medline import Shelf

    
def dbkeys(dbfile, keylist):
    """List the keys in a Berkeley database.
    
    @param dbfile: Path to Berkeley DB.
    @param keylist: Path to write database keys one per line.
    """
    d = db.DB()
    d.open(dbfile, None, db.DB_HASH, db.DB_RDONLY)
    cur = d.cursor()
    f = open(keylist, "w")
    rec = cur.first(dlen=0, doff=0)
    while rec is not None:
        f.write(rec[0]+"\n")
        rec = cur.next(dlen=0, doff=0)
    cur.close()
    d.close()
    f.close()
    

def pmid_dates(artdb, infile, outfile):
    """Get dates for PMIDs listed in L{infile}, writing PMID,date pairs
    to L{outfile} in increasing order of date.
    
    @param artdb: Path to Shelf with Article objects
    @param infile: Path to PubMed IDs (PMID lines)
    @param outfile: Path to write PMID YYYYMMDD lines to."""
    lines = []
    adb = Shelf.open(artdb, "r")
    input = open(infile, "r")
    for line in input:
        if line.startswith("#"): 
            continue
        pmid = int(line.split()[0])
        try:
            date = DateAsInteger(adb[str(pmid)].date_completed)
            lines.append((pmid,date))
        except KeyError, e:
            print e.message
    adb.close()
    input.close()
    lines.sort(key=lambda x:x[1])
    with open(outfile, "w") as f:
        for line in lines:
            f.write("%s %08d\n" % line)
    

def select_lines(infile, outfile, mindate="00000000", maxdate="99999999", N="0"):
    """Select random PMIDs from L{infile} and write them to L{outfile}.
    
    @param infile: Read PMID YYYYMMDD lines from this path.
    @param outfile: Write selected lines to this path.
    @param N: (string) Number of lines to output (N="0" outputs all matching)
    @param mindate, maxdate: Only consider YYYYMMDD strings between these
    @return: Selected lines as (PMID,YYYYMMDD) pairs of strings
    """
    lines = []
    N = int(N)
    input = open(infile, "r")
    for line in input:
        if line.startswith("#"): continue
        pmid, date = line.strip().split()
        if date >= mindate and date <= maxdate:
            lines.append((pmid,date))
    input.close()
    if N > 0:
        lines = random.sample(lines, N)
    lines.sort(key=lambda x:x[1])
    with open(outfile, "w") as f:
        for line in lines:
            f.write("%s %s\n" % line)
    return lines


if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger(logfile=False)
    locals()[sys.argv[1]](*sys.argv[2:])

