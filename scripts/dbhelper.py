#!/usr/bin/env python

"""Copy a FeatureDatabase to a FeatureStream

Usage::
    ./dbhelper.py listkeys <dbfile> <outfile>

List the keys in the Berkeley database <dbfile> writing them one per line to
<outfile>.

Usage::
    ./dbhelper restream <artdb> <featdb> <featstream>

Read Article dates from <artdb> (Shelf) and feature vectors from <featdb> 
(FeatureDatabase) to generate <featurestream> (FeatureStream).

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from __future__ import with_statement
from contextlib import closing
import sys
from bsddb import db
from mscanner.medline.FeatureDatabase import FeatureDatabase, FeatureStream
from mscanner.medline import Shelf

    
def list_db_keys(dbfile, outfile):
    """List the keys in infile, and print them one per line to outfile"""
    d = db.DB()
    d.open(dbfile, None, db.DB_HASH, db.DB_RDONLY)
    cur = d.cursor()
    f = open(outfile, "w")
    rec = cur.first(dlen=0, doff=0)
    while rec is not None:
        f.write(rec[0]+"\n")
        rec = cur.next(dlen=0, doff=0)
    cur.close()
    d.close()
    f.close()


def regenerate_stream(artdb, featdb, featstream):
    """Generate FeatureStream from the Shelf and FeatureDatabase"""
    adb = Shelf.open(artdb, "r")
    fd = FeatureDatabase(featdb, "r")
    fs = FeatureStream(open(featstream, "wb"))
    for pmid, art in adb.iteritems():
        fs.write(pmid, art.date_completed, fd[pmid])
    fs.close()
    fd.close()
    adb.close()


if __name__ == "__main__":
    action = sys.argv[1]
    if action == "listkeys":
        list_db_keys(sys.argv[2], sys.argv[3])
    elif action == "restream":
        regenerate_stream(sys.argv[2], sys.argv[3], sys.argv[4])

