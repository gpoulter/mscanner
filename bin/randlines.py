#!/usr/bin/env python

"""Selects a random subset of lines from one file and writes them to another

Usage::
    random-pmids.py [input file] [output file] [number of lines]

                                   

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import random
from path import path
import sys

def main(infile, outfile, nlines):
    if not infile.isfile():
        raise ValueError("%s does not exist" % infile)
    if outfile.isfile():
        raise ValueError("%s already exists!" % outfile)
    outfile.write_lines(random.sample(infile.lines(), nlines))

if __name__ == "__main__":
    main(path(sys.argv[1]), path(sys.argv[2]), int(sys.argv[3]))