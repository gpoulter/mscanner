#!env python

"""Convert feature database to feature stream...

Usage: ./db2stream.py <dbfile> <outputfile>

"""

from featuredb import FeatureDatabase, FeatureStream
import sys

if __name__ == "__main__":
    d = FeatureDatabase(sys.argv[1], 'r')
    s = FeatureStream(file(sys.argv[2], "wb"))
    for key, val in d.iteritems():
        s.write(key, val)
    d.close()
    s.close()
