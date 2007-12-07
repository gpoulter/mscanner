#!/usr/bin/env python

"""Calculates cross validation results for the MScanner paper.

Usage Example::
    python validate.py sample aidsbio pg07
    python validate.py iedb
"""

from path import path
import sys

from mscanner.configuration import rc
from mscanner.core import iofuncs
from mscanner.core.ValidationManager import SplitValidation, CrossValidation
from mscanner.medline.Databases import FeatureData, ArticleData

                                     
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
"aidsbio":       ("Paper/aidsbio.txt", "Paper/medline100k.txt"),
"pg07":          ("Paper/pg07_full.txt", "Paper/medline100k.txt"),
"radiology":     ("Paper/radiology.txt", "Paper/medline100k.txt"),
"control":       ("Paper/control.txt", "Paper/medline100k.txt"),
"gdsmall":       ("Test/gdsmall.txt", 30000),
}
"""Map data set to pair of (positive,negative) paths for PubMed IDs."""


def sample(*datasets):
    """Perform cross-validation analysis on the main data sets"""
    base = rc.working / "valid" / "Samples"
    if not base.exists(): base.mkdir()
    adata = ArticleData.Defaults()
    #fdata = FeatureData.Defaults_MeSH()
    fdata = FeatureData.Defaults_All()
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        pos, neg = dataset_map[dataset]
        if not isinstance(pos, path):
            pos = rc.corpora / pos
        if isinstance(neg, str):
            neg = rc.corpora / neg
        op = CrossValidation(base/dataset, dataset, adata, fdata)
        op.validation(pos, neg)
        op.report_validation()
        #op.report_predicted(1000, 10000, 16000000)


    
def trec():
    """Performs split-sample validation on the 2005 TREC Genomics track
    categorisation subtasks (Allele, Expression, GO, Tumor)."""
    base = rc.working / "valid" / "TREC"
    if not base.exists(): base.mkdir()
    adata = ArticleData.Defaults()
    fdata = FeatureData.Defaults_MeSH()
    #fdata = FeatureData.Defaults_All()
    ntest = rc.corpora / "TREC" / "NEG_test.txt"
    ntrain = rc.corpora / "TREC" / "NEG_train.txt"
    for ds, Ur in [("A",17.0), ("E",64.0), ("G",11.0), ("T",231.0)]:
        rc.utility_r = Ur
        dataset = "TREC_%s" % ds
        ptrain = rc.corpora / "TREC" / (ds + "train.txt")
        ptest = rc.corpora / "TREC" / (ds + "test.txt")
        op = SplitValidation(base/dataset, dataset, adata, fdata)
        op.validation(ptrain, ntrain, ptest, ntest)



def iedb():
    """Performs cross validation using the IEDB gold standard data set"""
    base = rc.working / "valid" / "IEDB_CV"
    if not base.exists(): base.mkdir()
    dataset = "IEDB CV MeSH DF4"
    rc.mincount = 4
    pos = rc.corpora / "IEDB" / "iedb-pos-dates.txt"
    neg = rc.corpora / "IEDB" / "iedb-neg-dates.txt"
    adata = ArticleData.Defaults()
    fdata = FeatureData.Defaults_MeSH()
    #fdata = FeatureData.Defaults_All()
    rc.utility_r = None
    op = CrossValidation(base/dataset, dataset, adata, fdata)
    op.validation(pos, neg)
    op.report_validation()
    
    
if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger(logfile=False)
    locals()[sys.argv[1]](*sys.argv[2:])


'''
def compare_score_methods():
    """Cross validation performance using the different feature score methods"""
    base = rc.working / "valid" / "ScoreMethod"
    if not base.exists(): base.mkdir()
    score_methods = ["scores_bayes", "scores_noabsence", "scores_rubin"]
    #train_P = rc.corpora / "pharmgkb-070205.txt"
    #train_N = rc.corpora / "medline07-100k.txt"
    train_P = rc.corpora / "genedrug-small.txt"
    train_N = rc.articlelist
    for method in score_methods:
        dataset = "pg07-" + method
        op = CrossValidation(base/dataset, dataset, env)
        op.validation(train_P, train_N)
        op.report_validation()
    
def issn_features():
    """Cross validation on for AIDSBio, excluding ISSN features. This shows that
    ISSN features produces better averaged precision in cross-validation."""
    rc.exclude_types = ["issn"]
    pos, neg = dataset_map["aidsbio"]
    op = CrossValidation(rc.working / "valid" / "aidsbio-noissn", 
                         "AIDSBio without ISSN")
    op.validation(rc.corpora / pos, rc.corpora / neg)
    op.report_validation()

def control_modality():
    """Cross validate Control against Medline100K, but without the features
    unique to Medline100K that cause the multi-modal score distributions."""
    pos, neg = dataset_map["control"]
    op = CrossValidation(rc.working / "valid" / "multimodal",
                         "Control Multi-modality")
    op.validation(rc.corpora / pos, rc.corpora / neg)
    op.report_validation()
'''
