"""Test program
"""

import numpy as n
from mscanner.scoring import iterCScores
from mscanner.featuredb import FeatureStream
import struct

if __name__ == "__main__":

    featscores = n.array([0.1, 5.0, 10.0, -5.0, -6.0])
    citations = [ (4, [4]), (4, [0,1,2]), (1,[0,2,3]), (2,[0,1]), (3,[1,2,3]) ]
    cite_fname = "citations.dat"

    # Write citations to disk
    fs = FeatureStream(file(cite_fname,"w"))
    for pmid, feats in citations:
        fs.write(pmid, n.array(feats, n.uint16))
    fs.close()

    # Calculate and print the scores
    for score, docid in iterCScores("cscore", cite_fname, len(citations), featscores, 5, None):
        print score, docid

