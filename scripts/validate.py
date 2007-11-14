#!/usr/bin/env python

"""Calculates cross validation results for the MScanner paper.

Choose the action by using a Python expression::
    python validate.py 'validate("aids-vs-500k")'
"""

from path import path
import sys

from mscanner.configuration import rc
from mscanner.core.ValidationManager import SplitValidation, CrossValidation
from mscanner.core import iofuncs
from mscanner.medline.Databases import Databases

                                     
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


dataset_map = {
"aidsbio":       ("Paper/aidsbio_2006.10.txt",    "Paper/medline100k_2007.02.txt"),
"pg07":          ("Paper/pharmgkb_2007.02.05.txt","Paper/medline100k_2007.02.txt"),
"radiology":     ("Paper/radiology_2007.02.txt",  "Paper/medline100k_2007.02.txt"),
"random10k":     ("Paper/random10k_2006.txt",     "Paper/medline100k_2007.02.txt"),
"gdsmall":       ("Test/gdsmall.txt",             100000),
}
"""Map data set to pair of (positive,negative) paths for PubMed IDs."""


def validate(*datasets):
    """Perform cross-validation analysis on the main data sets"""
    env = Databases()
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        pos, neg = dataset_map[dataset]
        if not isinstance(pos, path):
            pos = rc.corpora / pos
        if isinstance(neg, str):
            neg = rc.corpora / neg
        op = CrossValidation(rc.working / "valid" / dataset, dataset, env)
        op.validation(pos, neg)
        #op.report_validation()
        op.report_predicted(1000, 10000, 16000000)
    env.close()



def issn_features():
    """Perform cross validation on for AIDSBio vs 100k, excluding ISSN features.
    
    The results show that adding ISSN features produces better
    averaged precision in cross-validation."""
    rc.exclude_types = ["issn"]
    pos, neg = dataset_map["aidsbio"]
    op = CrossValidation(rc.working / "valid" / "aids-noissn", 
                         "AIDSBIO without ISSN features")
    op.validation(rc.corpora / pos, rc.corpora / neg)
    op.report_validation()
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
    rc.get_postmask = "mask_nonpositives"
    pos, neg = dataset_map["random10k"]
    op = CrossValidation(rc.working / "valid" / "random_modality",
                         "Removing bimodality in Random")
    op.validation(rc.corpora / pos, rc.corpora / neg)
    op.report_validation()
    op.env.close()


    
def trec2005_compare():
    """Performs split-sample validation on the trec GO subtask.
    
    The thing here is that TREC 2005 is an enriched background
    """
    env = Databases()
    ntest = rc.corpora / "TREC" / "NEG_test.txt"
    ntrain = rc.corpora / "TREC" / "NEG_train.txt"
    for ds, Ur in [("A",17.0), ("E",64.0), ("G",11.0), ("T",231.0)]:
        rc.utility_r = Ur
        dataset = "TREC_%s" % ds
        ptrain = rc.corpora / "TREC" / (ds + "train.txt")
        ptest = rc.corpora / "TREC" / (ds + "test.txt")
        op = SplitValidation(rc.working / "valid" / dataset, dataset, env)
        op.validation(ptrain, ntrain, ptest, ntest)
    env.close()



def wang2007_compare():
    """Performs cross validation using the data set from Wang2007"""
    env = Databases()
    rc.utility_r = None
    #for fbase in "ac", "allergen", "er", "other":
    for fbase in "combined",:
        pos = rc.corpora / "Wang2007" / ("%s_pos.txt" % fbase)
        neg = rc.corpora / "Wang2007" / ("%s_neg.txt" % fbase)
        dataset = "wang_%s" % fbase
        op = CrossValidation(rc.working / "valid" / dataset, dataset, env)
        op.validation(pos, neg)
        op.report_validation()
    env.close()
    
    

def compare_score_methods():
    """Compare cross validation performance using different feature
    score calculation methods on the PharmGKB data set."""
    score_methods = [
        "scores_offsetonly",
        "scores_withabsence",
        "scores_newpseudo",
        "scores_oldpseudo",
        "scores_rubin" ]
    env = Databases()
    #train_rel = rc.corpora / "pharmgkb-070205.txt"
    #train_irrel = rc.corpora / "medline07-100k.txt"
    train_rel = rc.corpora / "genedrug-small.txt"
    train_irrel = rc.articlelist
    for method in score_methods:
        dataset = "pg07-" + method
        op = CrossValidation(rc.working / "cmpscores" / dataset, dataset, env)
        op.validation(train_rel, train_irrel)
        op.report_validation()
    env.close()
    

if __name__ == "__main__":
    iofuncs.start_logger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression to execute"
    else:
        eval(sys.argv[1])
