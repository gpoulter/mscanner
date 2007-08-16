#!/usr/bin/env python

"""Small utility to compress shelved pickles in a Berkeley DB

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from mscanner.support import dbshelve

if __name__ == "__main__":
    print "OPENING"
    adb = dbshelve.open("articles_nocomp.db", "r", compress=False)
    cdb = dbshelve.open("articles.db", "w", compress=True)
    
    print "STARTING"
    next = 0
    for idx, (pmid, art) in enumerate(adb.iteritems()):
        if idx == next:
            print idx
            next += 1000
        cdb[pmid] = art
    print "ENDING"
    
    print "CLOSING"
    adb.close()
    cdb.close()
    print "DONE"
