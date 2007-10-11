#!/usr/bin/env python

"""Copy a FeatureDatabase to a FeatureStream

Usage::
    ./db2stream.py <dbfile> <outputfile>
    
This is to regenerate the FeatureStream (used by the cscore program to perform
fast queries), if it becomes corrupted but the FeatureDatabase is ok.

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import sys

def main(dbfile, streamfile):
    from mscanner.FeatureDatabase import FeatureDatabase, FeatureStream
    d = FeatureDatabase(dbfile, 'r')
    s = FeatureStream(open(streamfile, "wb"))
    for key, val in d.iteritems():
        s.write(key, val)
    d.close()
    s.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

