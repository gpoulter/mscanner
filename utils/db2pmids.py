#!env python

"""Extract a list of PubMed IDs from a FeatureDatabase

Usage: ./db2pmids.py <dbfile> <outputfile>

Iterate through <dbfile>, and print all the PubMed IDs therein to
<outputfile>.

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from medline import FeatureDatabase
import sys

if __name__ == "__main__":   
    d = FeatureDatabase(sys.argv[1], 'r')
    f = file(sys.argv[2], "w")
    for key in d:
        f.write("%d\n" % key)
    d.close()
    f.close()
