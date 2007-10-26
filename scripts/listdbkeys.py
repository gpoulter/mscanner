#!/usr/bin/env python

"""Extract a list of PubMed IDs from a FeatureDatabase

Usage::
    ./listdbkeys.py <dbfile> <outputfile>

Iterate through the Berkeley DB named <dbfile>, and print its keys one
per line to <outputfile>.

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import sys

def main(dbfile, outputfile):
    from bsddb import db
    d = db.DB() 
    d.open(dbfile, None, db.DB_HASH, db.DB_RDONLY)
    cur = d.cursor()
    f = open(outputfile, "w")
    rec = cur.first(dlen=0, doff=0)
    while rec is not None:
        f.write(rec[0]+"\n")
        rec = cur.next(dlen=0, doff=0)
    cur.close()
    d.close()
    f.close()

if __name__ == "__main__":   
    main(sys.argv[1], sys.argv[2])

