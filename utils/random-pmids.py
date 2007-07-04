#!python

"""Selects a random subset of lines from one file and writes them to another

Usage: random-pmids.py [input file] [output file] [number of lines]

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
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
