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
import codecs
from contextlib import closing
import logging
import numpy as nx
from path import path
import random
from pysqlite2 import dbapi2 as sqlite3
import sys

from mscanner.core import iofuncs
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
    

def pmid_dates(featdb, infile, outfile):
    """Get dates for PubMed IDs listed in L{infile}, writing PubMed ID and date
    pairs to L{outfile} in increasing order of date.
    @param featdb: Path to SQlite FeatureVectors.
    @param infile: Path to PubMed IDs (PMID lines).
    @param outfile: Path to write PMID YYYYMMDD lines to."""
    from mscanner.medline.FeatureVectors import FeatureVectors
    from mscanner.core.iofuncs import read_pmids
    fdb = FeatureVectors(featdb)
    pmids = read_pmids(infile)
    records = list((p,d) for (p,d,v) in fdb.get_records(pmids))
    fdb.close()
    records.sort(key=lambda x:x[1])
    with open(outfile, "w") as f:
        for pmid, date in records:
            f.write("%d %08d\n" % (pmid,date))
    

def select_lines(infile, outfile, mindate="00000000", maxdate="99999999", N="0"):
    """Select random PubMed IDs from L{infile} and write them to L{outfile}.
    
    @param infile: Read PubMed ID YYYYMMDD lines from this path.
    @param outfile: Write selected lines to this path.
    @param mindate, maxdate: Only consider YYYYMMDD strings between these values.
    @param N: (string) Number of lines to pint (N="0" outputs all matching lines)
    @return: Selected lines as (PMID,YYYYMMDD) pairs of strings
    """
    lines = []
    N = int(N)
    with open(infile, "r") as input:
        for line in input:
            if line.startswith("#"): continue
            pmid, date = line.strip().split()
            if date >= mindate and date <= maxdate:
                lines.append((pmid,date))
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

