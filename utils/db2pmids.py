#!env python

"""Extract PubMed IDs from the feature database

Usage: ./db2pmids.py <dbfile> <outputfile>

Iterate through <dbfile>, and print all the PubMed IDs therein to
<outputfile>.

"""

from medline import FeatureDatabase
import sys

d = FeatureDatabase(sys.argv[1], 'r')
f = file(sys.argv[1], "w")
for key in d:
    f.write("%d\n" % key)
d.close()
f.close()
