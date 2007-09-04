#!/usr/bin/env python

"""Convert a FeatureDatabase to a FeatureStream

Usage::
    ./db2stream.py <dbfile> <outputfile>

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want

"""

from mscanner import featuredb
import sys

def main():
    d = featuredb.FeatureDatabase(sys.argv[1], 'r')
    s = featuredb.FeatureStream(open(sys.argv[2], "wb"))
    for key, val in d.iteritems():
        s.write(key, val)
    d.close()
    s.close()

if __name__ == "__main__":
    main()