#!/usr/bin/env python

"""Creates exclude lists from the MeSH term tree

@author: Graham Poulter
                                   

Usage: make_excludes.py <file> [ <file> ]

The first file is the 'Mesh Trees' mentioned on
http://www.nlm.nih.gov/mesh/filelist.html, which can be downloaded. It
contains a list of MeSH terms with tree locations.

The second optional file is the 'ASCII Mesh Supplementary Records',
obtainable via 'MeSH in ASCII format' on the same page.  It contains
the supplementary MeSH terms.

Output 'mesh_exclude.cpickle', a set of all 'V' and 'Z' category terms
from the Mesh Trees, and also all main Supplementary Record headings.

Output 'mesh_synonyms.cpickle' mapping synonymous versions of
Supplementary Record headings to their main heading.

Place the output files in the pgmedline working directory

"""

import cPickle
import os
import re
import sys

def treeExcludes(excludes,tree):
    prefixes=["V","Z"]
    line=tree.readline()
    while line != "":
        line=line.strip()
        (term,code)=line.split(";")
        for p in prefixes:
            if code.startswith(p):
                excludes.add(term)
                break
        line=tree.readline()

def suppExcludes(synonyms,excludes,supp):
    line=supp.readline()
    while line != "":
        line=line.strip()
        if line=="*NEWRECORD":
            line=supp.readline()
            while line != "\n" and line != "":
                line=line.strip()
                if line.startswith("NM = "):
                    term=line[5:]
                    excludes.add(term)
                if line.startswith("SY = "):
                    syn=line[5:]
                    synonyms[syn]=term
                line=supp.readline()
        line=supp.readline()

if __name__=="__main__":
    if(len(sys.argv)<2):
        print __doc__
        sys.exit(1)
    if(not(os.path.exists(sys.argv[1]))):
        print __doc__
        sys.exit(1)
    excludes=set()
    treeExcludes(excludes,file(sys.argv[1],'r'))
    if(len(sys.argv)==3):
       if(not(os.path.exists(sys.argv[2]))):
           print __doc__
           sys.exit(1)
       synonyms=dict()
       suppExcludes(synonyms,excludes,file(sys.argv[2],'r'))
       #for (k,v) in synonyms.iteritems(): print k + " --> " + v
       cPickle.dump(synonyms,file("mesh_synonyms.cpickle",'wb'),protocol=2)
    #for t in excludes: print t
    cPickle.dump(excludes,file("mesh_excludes.cpickle",'wb'),protocol=2)
    print "Completed successfully"
