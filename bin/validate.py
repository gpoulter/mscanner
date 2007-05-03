#!/usr/bin/env python

"""Calculate performance statistics

validate.py dataset numnegs nfolds pseudocount alpha positives_path

                                   
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
from mscanner.utils import runMailer
from mscanner.validenv import ValidationEnvironment

def chooseValidation(datasets):
    env = ValidationEnvironment()
    rc.alpha = 0.5
    rc.cutoff = False
    rc.dodaniel = False
    rc.nfolds = 10
    rc.pseudocount = 0
    for dataset, pos, neg in [
        ("aids-vs-500k", "aids-bioethics-Oct06.txt", "medline07-500k.txt"),
        ("pg07-vs-500k", "pharmgkb-070205.txt", "medline07-500k.txt"),
        ("radiology-vs-500k", "daniel-radiology.txt", "medline07-500k.txt"),
        ("random10k-vs-500k", "random10k-06.txt", "medline07-500k.txt"), 
        ("gdsmall-vs-medline", "genedrug-small.txt", None) ]:
        if dataset not in datasets:
            continue
        rc.dataset = dataset
        pos = rc.corpora / pos
        neg = (rc.corpora / neg) if neg else rc.articlelist
        env.standardValidation(pos, neg)
        
def methodTests():
    env = ValidationEnvironment()
    rc.exclude_types = ["issn"]
    rc.dataset = "aids-noissn"
    pos = "aids-bioethics-Oct06.txt"
    neg = "medline07-500k.txt"
    env.standardValidation(rc.corpora / pos, rc.corpora / neg)

def compareDaniel():
    pg04 = rc.corpora / "pharmgkb-2004.txt"
    m30k = rc.corpora / "medline07-30k.txt"
    m500k = rc.corpora / "medline07-500k.txt"
    pg07 = rc.corpora / "pharmgkb-070205.txt"
    env = ValidationEnvironment()
    rc.alpha = 0.5
    rc.pseudocount = 0
    rc.nfolds = 10

    rc.dataset = "pg04-vs-30k"
    rc.dodaniel = False
    rc.exclude_types = None
    env.standardValidation(pg04, m30k)

    rc.dataset = "pg04-vs-30k-dan"
    rc.exclude_types = ["issn"]
    rc.dodaniel = True
    env.standardValidation(pg04, m30k)
    
    rc.dataset == "pg04-vs-500k"
    rc.dodaniel = False
    rc.exclude_types = None
    env.standardValidation(pg04, m500k)
    
    rc.dataset = "pg07-vs-30k"
    rc.dodaniel = False
    rc.exclude_types = None
    env.standardValidation(pg07, m30k)

    rc.dataset == "pg07-vs-500k-dan"
    rc.dodaniel = True
    rc.exclude_types = ["issn"]
    env.standardValidation(pg07, m30k)
    
    rc.dataset == "pg07-vs-500k-noissn"
    rc.dodaniel = False
    rc.exclude_types = ["issn"]
    env.standardValidation(pg07, m500k)
        
    #elif dataset == "mscanner-vs-500k":
    #    pos = "mscanner-bibliography.txt"
    #    neg = "medline07-500k.txt"
    #elif dataset == "pg07-vs-med07":
    #    pos = "pharmgkb-070205.txt"
    #    neg = articlelist
    #elif dataset == "gdsmall-vs-sample":
    #    pos = "genedrug-small.txt"
    #    neg = c.articlelist

def scriptmain(*args):
    """Meant to be called with *sys.argv[1:]"""
    if len(args) == 0:
        raise ValueError("Please give dataset code or CGI parameters")
    elif len(args) == 1:
        chooseValidation(set([args[0]]))
    elif len(args) > 1:
        rc.dataset = sys.argv[0]
        rc.numnegs = int(sys.argv[1])
        rc.nfolds = int(sys.argv[2])
        rc.pseudocount = float(sys.argv[3])
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
