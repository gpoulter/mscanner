#!/usr/bin/env python

"""Calculate performance statistics

Usage as backend for CGI:

  validate.py dataset numnegs nfolds alpha positives_path

Usage for experiments (executable python string):

  validate.py 'paperTests("aids-vs-500k")'

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

from path import path
import sys

from mscanner.configuration import rc, initLogger
from mscanner.statusfile import runMailer
from mscanner.validenv import ValidationEnvironment

dataset_map = {
    "aids-vs-100k": ("aids-bioethics-Oct06.txt", "medline07-100k.txt"),
    "pg07-vs-100k": ("pharmgkb-070205.txt", "medline07-100k.txt"),
    "radiology-vs-100k": ("daniel-radiology.txt", "medline07-100k.txt"),
    "random10k-vs-100k": ("random10k-06.txt", "medline07-100k.txt"), 
    "gdsmall-vs-sample": ("genedrug-small.txt", rc.articlelist) ,
    }

def paperTests(*datasets):
    """Cross-validation tests for the publication"""
    env = ValidationEnvironment()
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        rc.dataset = dataset
        pos, neg = dataset_map[dataset]
        if not isinstance(pos, path):
            pos = rc.corpora / pos
        if not isinstance(neg, path):
            neg = rc.corpora / neg
        env.standardValidation(pos, neg)
        
def issnTest():
    """Output cross-validatin results on AIDSBio when ISSN features are 
    excluded.  This will show how much including ISSN helps. 
    """
    env = ValidationEnvironment()
    rc.exclude_types = ["issn"]
    rc.dataset = "aids-noissn"
    pos = "aids-bioethics-Oct06.txt"
    neg = "medline07-500k.txt"
    env.standardValidation(rc.corpora / pos, rc.corpora / neg)

def compareDaniel():
    """Make comparisons against the Rubin2005 method.  
    
    @note: We test Daniel's smoothing versus Bayes prior for old and
    new pharmacogenetics data, against 30k or 500k negatives.
    
    """
    pg04 = rc.corpora / "pharmgkb-2004.txt"
    m30k = rc.corpora / "medline07-30k.txt"
    m500k = rc.corpora / "medline07-500k.txt"
    pg07 = rc.corpora / "pharmgkb-070205.txt"
    env = ValidationEnvironment()
    
    # Compare on PG04
    rc.dataset = "pg04-vs-30k"
    rc.exclude_types = None
    rc.getFrequencies = "getProbabilitiesBayes"
    env.standardValidation(pg04, m30k)

    rc.dataset = "pg04-vs-30k-dan"
    rc.exclude_types = ["issn"]
    rc.getFrequencies = "getProbabilitiesRubin"
    env.standardValidation(pg04, m30k)
    
    # Compare daniel vs default on PG07
    rc.dataset = "pg07-vs-30k"
    rc.exclude_types = None
    rc.getFrequencies = "getProbabilitiesBayes"
    env.standardValidation(pg07, m30k)

    rc.dataset == "pg07-vs-500k-dan"
    rc.exclude_types = ["issn"]
    rc.getFrequencies = "getProbabilitiesRubin"
    env.standardValidation(pg07, m30k)
    
    rc.dataset == "pg07-vs-500k-noissn"
    rc.exclude_types = ["issn"]
    rc.getFrequencies = "getProbabilitiesBayes"
    env.standardValidation(pg07, m500k)
    
def random_bimodality():
    """Normally Random10K under cross validation against other random
    citations has multiple modes in the article scores.
    
    This tests whether eliminating features that occur only once or
    twice in the data is effective at removing the bimodality (since
    such features are hypothesised to cause a spike of negative
    feature scores in training, and random occurrences in testing
    citations causes their scores to be shifted left from the mean).    
    """
    rc.dataset = "random10k-vs-100k-mod"
    rc.nfolds = 10
    rc.getPostMask = "maskLessThanThree"
    env = ValidationEnvironment()
    pos = "random10k.txt"
    neg = "medline07-100k.txt"
    #pos = "genedrug-small.txt"
    #neg = rc.articlelist
    env.standardValidation(rc.corpora / pos, rc.corpora / neg)

def scriptmain(*args):
    """Meant to be called with *sys.argv[1:], such that 
    the positional arguments of this function are the arguments to
    the program."""
    if len(args) == 0:
        raise ValueError("Please give dataset code or CGI parameters")
    elif len(args) == 1:
        eval(args[0])
    elif len(args) > 1:
        rc.dataset = sys.argv[0]
        rc.numnegs = int(sys.argv[1])
        rc.nfolds = int(sys.argv[2])
        rc.alpha = float(sys.argv[4])
        positives_path = path(sys.argv[5])
        negatives_path = rc.valid_report_dir / rc.report_negatives
        negatives_path.write_lines(
            random.sample(rc.articlelist.lines(), rc.numnegs))
        try:
            env = ValidationEnvironment()
            env.standardValidation(positives_path, negatives_path)
        finally:
            runMailer(rc.smtpserver, rc.emails_path)
    
if __name__ == "__main__":
    initLogger()
    scriptmain(*sys.argv[1:])
