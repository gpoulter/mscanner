#!env python

"""Filter positive and negative corpora

@author: Graham Poulter
                                   

Usage: fixcopora.py <positive> <negative> <bad>

Command line arguments specify files listing inputs of positive,
negative and bad PubMed IDs.  

<bad> specifies a file of 'broken' PMIDs to be removed from <positive>
and <negative> files.

Output is lists of sorted, unique positive and negative articles in
pos.txt and neg.txt.  Articles found <bad> list are excluded from
pos.txt, and negative articles found in <positive> are left out of the
neg.txt.

"""

import sys, os.path

def read_positives(fname):
    lines=file(fname,'r').readlines()
    result=set()
    for l in lines:
        terms=l.strip().split()
        result.add(terms[0])
    return result

def read_pmids(fname):
    lines = file(fname,'r').readlines()
    return set([l.strip() for l in lines])

if __name__=="__main__":

    if len(sys.argv) != 4:
        raise ValueError("Incorrect number of command-line arguments.")

    if not os.path.exists(sys.argv[1]) or \
       not os.path.exists(sys.argv[2]) or \
       not os.path.exists(sys.argv[3]):
        raise ValueError("Could not locate one of the specified files")

    pos_in=read_positives(sys.argv[1])
    neg_in=read_pmids(sys.argv[2])
    bad=read_pmids(sys.argv[3])

    pos_out=[p for p in pos_in if p not in bad]
    pos_out.sort(key=int)
    
    neg_out=[p for p in neg_in if p not in pos_in|bad]
    neg_out.sort(key=int)

    file("pos.txt",'w').writelines(p+"\n" for p in pos_out);
    file("neg.txt",'w').writelines(p+"\n" for p in neg_out);
