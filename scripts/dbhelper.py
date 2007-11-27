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
from mscanner.medline.FeatureStream import FeatureStream, Date2Integer
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


def featuredata():
    """Regenerate the feauture databases from the article database. Please
    move the old streams out of the way first!"""
    updater = Updater.Defaults()
    print "Updating MeSH-features data"
    updater.fdata_mesh.add_articles((str(a.pmid), 
        Date2Integer(a.date_completed), 
        a.mesh_features())
        for a in updater.adata.artdb.itervalues())
    print "Updating all-features data"
    updater.fdata_all.add_articles((str(a.pmid), 
        Date2Integer(a.date_completed), 
        a.all_features(updater.stopwords))
        for a in updater.adata.artdb.itervalues())

    

def featurestream(artdb, featdb, featstream, bits):
    """Regenerate a FeatureStream from an article db and FeatureDatabase.
    @param artdb: Path to Shelf with Article objects
    @param featdb: Path to FeatureDatabase
    @param featstream: Path to write re-generated FeatureStream to
    @param bits: Either "16" or "32" for number of bits in a feature
    """
    ftype = {"16":nx.uint16,"32":nx.uint32}[bits]
    adb = Shelf.open(artdb, "r")
    fd = FeatureDatabase(featdb, "r", ftype=ftype)
    if path(featstream).isfile():
        path(featstream).remove()
    fs = FeatureStream(featstream, "ab", ftype=ftype)
    for idx, (pmid, art) in enumerate(adb.iteritems()):
        if idx % 10000 == 0:
            print "Completed %d" % idx
        fs.additem(pmid, art.date_completed, fd[pmid])
    fs.close()
    fd.close()
    adb.close()


def articlelist(artdb, artlist):
    """Regenerate the list of articles from the article database
    @param artdb: Path to Shelf with Article objects
    @param articlelist: Path to write PubMed IDs and YYYYMMDD lines to."""
    adb = Shelf.open(artdb, "r")
    f = open(artlist, "w")
    for idx, (pmid, art) in enumerate(adb.iteritems()):
        if idx % 10000 == 0:
            print "Completed %d" % idx
        f.write("%s %08d\n" % (pmid, Date2Integer(art.date_completed)))
    f.close()
    adb.close()
    
    
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
            date = Date2Integer(adb[str(pmid)].date_completed)
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

