#!/usr/bin/env python

"""Tests different scoring methods to see whether there is a
substantial difference in performance"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import sys

from mscanner.configuration import rc, start_logger
from mscanner.medline.Databases import Databases
from mscanner.core.Validator import CrossValidator


score_methods = [
    "scores_offsetonly",
    "scores_withabsence",
    "scores_newpseudo",
    "scores_oldpseudo",
    "scores_rubin" ]


def do_comparisons():
    env = Databases()
    #train_rel = rc.corpora / "pharmgkb-070205.txt"
    #train_irrel = rc.corpora / "medline07-100k.txt"
    train_rel = rc.corpora / "genedrug-small.txt"
    train_irrel = rc.articlelist
    for method in score_methods:
        rc.dataset = "pg07-" + method
        v = CrossValidation(rc.working / "cmpscores" / rc.dataset, env)
        v.validation(train_rel, train_irrel)
    env.close()
    

if __name__ == "__main__":
    start_logger()
    do_comparisons()
