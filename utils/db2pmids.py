#!/usr/bin/env python

"""Extract PubMed IDs from the feature database

Usage: ./extractpmids.py <outputfile>

Iterate through the configured feature database, and print all the
PubMed IDs therein to <outputfile>.

"""

import configuration as c
from medline import FeatureDatabase
import sys

d = FeatureDatabase(c.medline.featuredb, 'r')
f = file(sys.argv[1], "w")
for key in d:
    f.write("%d\n" % key)
d.close()
f.close()
