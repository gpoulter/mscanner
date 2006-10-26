#!env python

"""Extract PubMed IDs from the feature database

Usage: ./db2pmids.py <pmidfile> <featdbfile> <outputfile>

Iterate through <dbfile>, and print all the PubMed IDs therein to
<outputfile>.

"""

import sys
import dbshelve
from medline import FeatureDatabase

pmids = [ int(line) for line in file(sys.argv[1]) if not line.startswith("#") ]
fdb = FeatureDatabase(sys.argv[2], 'rw')
adb = dbshelve.open(sys.argv[3], 'rw')
ofile = file("narticles.txt","w")

for pmid in pmids:
    if str(pmid) not in adb:
        print "%d not in articles" % pmid
        if pmid in fdb:
            print "deleting %d from features" % pmid
            fdb.delitem(pmid)
    if pmid not in fdb:
        print "%d not in features" % pmid
        if str(pmid) in adb:
            print "deleting %d from articles" % pmid
            adb.delitem(pmid)
    if str(pmid) in adb and pmid in fdb:
        ofile.write("%d\n"%pmid)

fdb.close()
adb.close()
ofile.close()
