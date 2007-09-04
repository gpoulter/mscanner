#!/usr/bin/env python

"""Extract PubMed IDs from GeneOntology annotation files

Usage::
    go2pmids.sh <file1.gz> <file2.gz> ...

Files should be .gz (gzip) files containing GeneOntology annotations.
The program prints a sorted list of unique PubMed IDs (some of which 
are Medline UIs or gibberish) to standard output.

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import sys
import os
import re
import gzip

def parseGeneOntology(fname,pmids):
    pmid_finder = re.compile(r'PMID:([0-9]{1,8})')
    input = gzip.open(fname)
    for line in input:
        result = pmid_finder.search(line)
        if result is not None:
            pmids.add(result.group(1))
    input.close()
            
def main():
    if(len(sys.argv)==1):
        print __doc__
        sys.exit(0)
    files = sys.argv[1:]
    for f in files:
        if not os.path.exists(f):
            print "Error: File %s does not exist" % (f,)
            sys.exit(1)
    pmids = set()
    for f in files:
        parseGeneOntology(f,pmids)
    pmids = list(pmids)
    pmids.sort(key=int)
    for p in pmids:
        print p

if __name__=="__main__":
    main()