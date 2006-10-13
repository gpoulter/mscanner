#!/usr/bin/env python

"""Extracts PubMed IDs from GeneOntology annotation files

@author: Graham Poulter
                                   

Usage: parse_geneontology.sh <file1.gz> <file2.gz> ...

Files should be .gz (gzip) files containing GeneOntology annotations.
The program prints a sorted list of unique PubMed IDs from the to
standard output.

"""

import sys, os, re, gzip

def parseGeneOntology(fname,pmids):
    pmid_finder=re.compile(r'PMID:([0-9]{1,8})')
    input=gzip.open(fname)
    for line in input:
        result=pmid_finder.search(line)
        if result is not None:
            pmids.add(result.group(1))

if __name__=="__main__":
    if(len(sys.argv)==1):
        print __doc__
        sys.exit(0)
    files=sys.argv[1:]
    for f in files:
        if not os.path.exists(f):
            print "Error: File %s does not exist" % (f,)
            sys.exit(1)
    pmids=set()
    for f in files:
        parseGeneOntology(f,pmids)
    pmids=list(pmids)
    pmids.sort(key=int)
    for p in pmids:
        print p
