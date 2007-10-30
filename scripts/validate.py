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
from mscanner.medline.Databases import Databases
from mscanner.core.ValidationManager import SplitValidation, ValidationManager


dataset_map = {
"aidsbio":       ("Paper/aidsbio_2006.10.txt",    "Paper/medline100k_2007.02.txt"),
"pg07":          ("Paper/pharmgkb_2007.02.05.txt","Paper/medline100k_2007.02.txt"),
"radiology":     ("Paper/radiology_2007.02.txt",  "Paper/medline100k_2007.02.txt"),
"random10k":     ("Paper/random10k_2006.txt",     "Paper/medline100k_2007.02.txt"),
"gdsmall":       ("Test/gdsmall.txt",             None),
"trec_go_bg" :   ("TREC/trec_go_pos_both.txt",    "TREC/trec_go_bg_both.txt"),
"trec_go_100k" : ("TREC/trec_go_pos_both.txt",    "Paper/medline100k_2007.02.txt"),
}
"""Map data set to pair of (positive,negative) paths for PubMed IDs."""


def validate(*datasets):
    """Perform cross-validation analysis for the Mscanner publication"""
    env = Databases()
    rc.numnegs = 100000
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        rc.dataset = dataset
        pos, neg = dataset_map[dataset]
        if not isinstance(pos, path):
            pos = rc.corpora / pos
        if neg is not None and not isinstance(neg, path):
            neg = rc.corpora / neg
        op = ValidationManager(rc.working / "valid" / rc.dataset, env)
        op.validation(pos, neg, rc.nfolds)
    env.close()


def issn_validation():
    """Perform cross validation on for AIDSBio vs 100k, excluding ISSN features.
    
    The results show that keeping the ISSN features produces better
    cross-validation performance."""
    rc.exclude_types = ["issn"]
    rc.dataset = "aids-noissn"
    pos, neg = dataset_map["aidsbio"]
    op = ValidationManager(rc.working / "valid" / rc.dataset)
    op.validation(rc.corpora / pos, rc.corpora / neg, rc.nfolds)
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
    rc.dataset = "random_modality"
    rc.nfolds = 10
    rc.get_postmask = "mask_nonpositives"
    pos, neg = dataset_map["random10k"]
    op = ValidationManager(rc.working / "valid" / rc.dataset)
    op.validation(rc.corpora / pos, rc.corpora / neg, rc.nfolds)
    op.env.close()
    
    
def trec2005_compare():
    """Performs split-sample validation on the trec GO subtask"""
    env = Databases()
    rc.nfolds = 1
    ntest = rc.corpora / "TREC" / "NEG_test.txt"
    ntrain = rc.corpora / "TREC" / "NEG_train.txt"
    for ds, Ur in [("A",17.0), ("E",64.0), ("G",11.0), ("T",231.0)]:
        rc.dataset = "TREC_%s" % ds
        rc.utility_r = Ur
        ptrain = rc.corpora / "TREC" / (ds + "train.txt")
        ptest = rc.corpora / "TREC" / (ds + "test.txt")
        op = SplitValidation(rc.working / "valid" / rc.dataset, env)
        op.validation(ptrain, ntrain, ptest, ntest)
    env.close()


def wang2007_compare():
    env = Databases()
    for fbase in "ac", "allergen", "er", "other":
        pos = rc.corpora / "Wang2007" / ("wang_%s_pos.txt" % fbase)
        neg = rc.corpora / "Wang2007" / ("wang_%s_neg.txt" % fbase)
        rc.dataset = "wang_%s" % fbase
        rc.nfolds = 10
        rc.utility_r = None
        op = ValidationManager(rc.working / "valid" / rc.dataset, env)
        op.validation(pos, neg, rc.nfolds)
    env.close()
    


if __name__ == "__main__":
    start_logger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression"
    else:
        eval(sys.argv[1])
