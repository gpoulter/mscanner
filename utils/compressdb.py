#!env python

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
