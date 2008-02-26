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
import sqlite3
import sys

from mscanner.core import iofuncs
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
    from mscanner.medline.FeatureStream import DateAsInteger
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


def upgrade_featmap(infile, outfile):
    """Read in the old FeatureMapping text format and write the SQLite version"""
    with closing(codecs.open(infile, "rb", "utf-8")) as input:
        input.readline() # Get past article count
        with closing(sqlite3.connect(outfile)) as con:
            con.execute("""CREATE TABLE IF NOT EXISTS fmap (
              id INTEGER PRIMARY KEY,
              type TEXT, name TEXT, count INTEGER,
              UNIQUE(type,name) )""")
            for fid, line in enumerate(input):
                fname, ftype, count = line.split("\t")
                count = int(count)
                con.execute("INSERT INTO fmap VALUES(?,?,?,?)", 
                            (fid, ftype, fname, count))
            con.commit()


def upgrade_featdb(infile, outfile):
    """Read in FeatureStream write SQLite version"""
    from mscanner.medline.FeatureDatabase import FeatureVectors
    from mscanner.medline.FeatureStream import FeatureStream
    with closing(FeatureStream(path(infile),rdonly=True)) as fs:
        fv = FeatureVectors(outfile)
        count = 0
        for pmid, date, vector in fs.iteritems():
            if pmid in fv:
                logging.warning("PMID %d was found already!\n")
            else:
                fv.add_record(pmid, date, vector)
                #assert fv.get_vector(pmid) == vector
            count += 1
            if count % 10000 == 0:
                logging.debug("Written %d articles.", count)
                fv.con.commit()
        fv.con.commit()
        logging.debug("Wrote %d (%d) articles in total.", count, len(fv))


if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger(logfile=False)
    locals()[sys.argv[1]](*sys.argv[2:])

