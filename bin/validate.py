#!/usr/bin/env python

"""Calculates cross validation results for the MScanner paper.

Choose the action by using a Python expression::
    python validate.py 'validate("aids-vs-500k")'
"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


from path import path
import sys

from mscanner.configuration import rc, start_logger
from mscanner import scorefile, validenv


dataset_map = {
    "aids-vs-100k": ("aids-bioethics-Oct06.txt", "medline07-100k.txt"),
    "pg07-vs-100k": ("pharmgkb-070205.txt", "medline07-100k.txt"),
    "radiology-vs-100k": ("daniel-radiology.txt", "medline07-100k.txt"),
    "random10k-vs-100k": ("random10k-06.txt", "medline07-100k.txt"), 
    "gdsmall-vs-sample": ("genedrug-small.txt", rc.articlelist) ,
    "test-vs-sample": ("testing-random.txt", rc.articlelist),
}
"""Map data set to pair of (positive,negative) paths for PubMed IDs."""


def validate(*datasets):
    """Perform cross-validation analysis for the Mscanner publication"""
    env = scorefile.Databases()
    rc.numnegs = 1000
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        rc.dataset = dataset
        pos, neg = dataset_map[dataset]
        if not isinstance(pos, path):
            pos = rc.corpora / pos
        if neg is not None and not isinstance(neg, path):
            neg = rc.corpora / neg
        op = validenv.Validation(rc.working / "valid" / rc.dataset, env)
        op.validation(pos, neg)
    env.close()


def issn_validation():
    """Perform cross validation on for AIDSBio vs 100k, excluding ISSN features.
    
    The results show that keeping the ISSN features produces better
    cross-validation performance."""
    rc.exclude_types = ["issn"]
    rc.dataset = "aids-noissn"
    pos = "aids-bioethics-Oct06.txt"
    neg = "medline07-100k.txt"
    op = validenv.Validation(rc.working / "valid" / rc.dataset)
    op.validation(rc.corpora / pos, rc.corpora / neg)
    op.env.close()


def random_bimodality():
    """Cross validation of Random10K against Medline100K results
    in multiple modes in the article score distributions.
    
    The distribution becomes unimodal about zero if we eliminate features that
    do not occur in positive citations. There are thousands of such features
    all of which have a negative feature scores of roughly the same value (-9).
    Random occurrences of 0,1,2 or 3 of those -9 features in test citations
    causes the score of that citation to be shifted by -9, -18, -27 from
    zero."""
    rc.dataset = "random10k-vs-100k-mod"
    rc.nfolds = 10
    rc.post_masker = "maskNonPositives"
    pos = "random10k-06.txt"
    neg = "medline07-100k.txt"
    op = validenv.Validation(rc.working / "valid" / rc.dataset)
    op.validation(rc.corpora / pos, rc.corpora / neg)
    op.env.close()


if __name__ == "__main__":
    start_logger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression"
    else:
        eval(sys.argv[1])
