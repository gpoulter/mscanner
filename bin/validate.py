#!/usr/bin/env python

"""Cross-validation performance analysis

Choose operation by executing a Python expression::
    python validate.py 'paperTests("aids-vs-500k")'
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

from mscanner.configuration import rc, initLogger
from mscanner.validenv import ValidationEnvironment


def paperTests(*datasets):
    """Cross-validation tests for the publication"""
    dataset_map = {
    "aids-vs-100k": ("aids-bioethics-Oct06.txt", "medline07-100k.txt"),
    "pg07-vs-100k": ("pharmgkb-070205.txt", "medline07-100k.txt"),
    "radiology-vs-100k": ("daniel-radiology.txt", "medline07-100k.txt"),
    "random10k-vs-100k": ("random10k-06.txt", "medline07-100k.txt"), 
    "gdsmall-vs-sample": ("genedrug-small.txt", rc.articlelist) ,
    "test-vs-sample": ("testing-random.txt", rc.articlelist),
    }
    env = ValidationEnvironment()
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
        env.standardValidation(pos, neg)


def issnTest():
    """Output results on AIDSBio when ISSN features are excluded. Indicates
    the degree of performance improvement due to ISSNs."""
    env = ValidationEnvironment()
    rc.exclude_types = ["issn"]
    rc.dataset = "aids-noissn"
    pos = "aids-bioethics-Oct06.txt"
    neg = "medline07-500k.txt"
    env.standardValidation(rc.corpora / pos, rc.corpora / neg)


def random_bimodality():
    """Normally Random10K under cross validation against other random
    citations has multiple modes in the article scores.
    
    This tests whether eliminating features that do not occur in positive
    citations removes the bimodality (since such features cause a spike of
    negative feature scores in training, and random occurrences in testing
    citations causes their scores to be shifted left from the mean)."""
    rc.dataset = "random10k-vs-100k-mod"
    rc.nfolds = 10
    rc.post_masker = "maskNonPositives"
    env = ValidationEnvironment()
    pos = "random10k-06.txt"
    neg = "medline07-100k.txt"
    env.standardValidation(rc.corpora / pos, rc.corpora / neg)


if __name__ == "__main__":
    initLogger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression"
    else:
        eval(sys.argv[1])
