#!python

"""Selects a random subset of lines from one file and writes them to another

Usage: random-pmids.py [input file] [output file] [number of lines]
"""

import random
from path import path
import sys

if __name__ == "__main__":

    infile = path(sys.argv[1])
    
    if not infile.isfile():
        raise ValueError("%s does not exist" % infile)
    
    outfile = path(sys.argv[2])
    
    if outfile.isfile():
        raise ValueError("%s already exists!" % outfile)
    
    nlines = int(path(sys.argv[3]))
    
    outfile.write_lines(random.sample(infile.lines(), nlines))
