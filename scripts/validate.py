#!/usr/bin/env python

"""Calculates cross validation results for the MScanner paper.

Usage Example::
    python validate.py sample aidsbio pg07
    python validate.py iedb
"""

import logging
import numpy as nx
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


# Top-level directory for all outputs
base = rc.root / "results"

# Data set codes to input paths
dataset_map = {
    "aidsbio":       ("Paper/aidsbio.txt", "Paper/medline100k.txt"),
    "pg07":          ("Paper/pg07_full.txt", "Paper/medline100k.txt"),
    "radiology":     ("Paper/radiology.txt", "Paper/medline100k.txt"),
    "control":       ("Paper/control.txt", "Paper/medline100k.txt"),
    "gdsmall":       ("Test/gdsmall.txt", 30000),
    "iedb":          ("IEDB/iedb-pos-dates.txt", "IEDB/iedb-neg-dates.txt"),
}

# Map feature space codes to full spaces and excluded classes
spaces = {
    "iedbword":     ("feats_iedb_word", []),
    "iedbconcat":   ("feats_iedb_concat", []),
    
    "wmqia":        ("feats_wmqia", []),
    "wmqiafilt":    ("feats_wmqia_filt", []),
    "word":         ("feats_wmqia", ["mesh","qual","issn","a"]),
    "author":       ("feats_wmqia", ["mesh","qual","issn","w"]),
    
    "meshqi":       ("feats_mesh_qual_issn", []),
    "meshq":        ("feats_mesh_qual_issn", ["issn"]),
    "meshnq":       ("feats_mesh_qual_issn", ["qual","issn"]),
    "qual":         ("feats_mesh_qual_issn", ["mesh","issn"]),
    "issn":         ("feats_mesh_qual_issn", ["mesh","qual"]),
}

# Cache of loaded feature space databases (featurespace->FeatureData)
fdata_cache = {}

# Cache of loaded data sets (path->array of PMIDs)
pmids_cache = {}

class ResultsTable:
    """Class for tabulating runtime parameters and performance statistics
    results.
    
    @ivar fname: Name of file to which CSV table is written.
    @ivar stream: File object that its written to.
    """
    
    def __init__(self, fname, append=False):
        self.fname = fname
        if not fname.parent.exists():
            fname.parent.makedirs()
        write_title = not (fname.exists() and append)
        self.stream = open(fname, "a" if append else "w")
        if write_title:
            self.stream.write("Name, ROC, ROC SD, AvPrec, BEP, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\n")
        
    def close(self):
        self.stream.close()
        
    def add_results(self, vmanager):
        """Add a row of results from a L{ValidationManager}."""
        v = vmanager.metric_vectors
        results = [vmanager.dataset, v.W, v.W_stderr, v.AvPrec, v.breakeven]
        fmt = "%s, %.4f, %.4f, %.3f, %.3f" + (", %.2f"*11) + "\n"
        self.stream.write(fmt % tuple(results + v.precision_11))
        self.stream.flush()


def get_dataset(dataset, fdata):
    """Given a data set name, return (pos, neg) pair of arrays of PMIDs for
    positive and negative articles."""
    if dataset not in dataset_map:
        raise ValueError("Invalid Data Set %s" % dataset)
    # Get pair of incomplete paths
    dspaths = list(dataset_map[dataset])
    results = []
    for dspath in dspaths:
        if isinstance(dspath, str):
            # Complete the path, and cache loaded PMIDs for later
            try:
                pmids = pmids_cache[dspath]
            except KeyError:
                pmids, notfound, exclude = \
                     iofuncs.read_pmids_careful(rc.corpora / dspath, fdata.featuredb)
                pmids_cache[dspath] = pmids
            results.append(pmids)
        else:
            # Let ValidationManager load the data
            results.append(dspath)
    # Convert list back to array
    return tuple(results)


def base_valid(outdir, dataset, featurespace, tab=None, df=None, ig=None):
    """General purpose cross-validation.

    @param outdir: Directory in which to place output dir (relative path from
    base is used to derive the title of the dataset)

    @param dataset: Input corpus to use
    @param featurespace: Name of feature space to use
    @param df: Minimum Document Frequency (rc.mincount)
    @param ig: Minimum Information Gain (rc.min_infogain)
    @param tab: Instance of L{ResultsTable}
    """
    # Skip the task if it's already been done
    #if outdir.exists(): return
    # Make output directory and calculate title
    if not outdir.parent.exists(): 
        outdir.parent.makedirs()
    # Create the data set title
    title = "_".join(base.relpathto(outdir).splitall()[1:])
    # Set document frequency cutoff
    if df is not None:
        rc.mincount = df
    # Set information gain cutoff
    if ig is not None:
        rc.min_infogain = ig
    # Choose feature types
    if featurespace == "feats_mesh_qual_issn":
        ftype = nx.uint16
    else:
        ftype = nx.uint32
    # Cache loaded feature spaces (saves a few seconds)
    if featurespace not in fdata_cache:
        fdata_cache[featurespace] = FeatureData.Defaults(featurespace, ftype)
    fdata = fdata_cache[featurespace]
    # Perform cross validation
    op = CrossValidation(outdir, title, adata, fdata)
    op.validation(*get_dataset(dataset, fdata))
    op.report_validation()
    if tab is not None:
        tab.add_results(op)


def compare_prior(ds):
    """Different priors"""
    rc.mincount = 1
    rc.min_infogain = 0
    rc.positives_only = False
    s = base / "PRIOR" / "df1_ig0_all" / ds / ("sd%d"%rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for fs in ["meshq","word","wmqia"]:
        fspace, rc.class_mask = spaces[fs]
        for method in ["bgfreq", "laplace_split", "laplace"]:
            rc.scoremethod = "scores_" + method
            base_valid(s/fs/method, ds, fspace, tab)


def compare_featselection(ds):
    """Different feature selection approaches."""
    rc.scoremethod = "scores_laplace_split"
    s = base / "SELN" / "lsplit" / ds / ("sd%d"%rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for fs in ["meshq","word","wmqia"]:
        fspace, rc.class_mask = spaces[fs]
        # Information Gain
        for ig in [0, 1e-5, 2e-5, 1e-4]:
            base_valid(s/fs/("ig%.1f_df1_all"%(ig*1e5)), ds, fspace, tab, df=1, ig=ig)
        # Document Frequency
        for df in [0, 4, 8]:
            base_valid(s/fs/("df%d_ig0_all"%df), ds, fspace, tab, df=df, ig=0)
        # Only if present in a relevant example
        rc.positives_only = True
        base_valid(s/fs/"pos_df1_ig0", ds, fspace, tab, df=1, ig=0)
        rc.positives_only = False


def compare_featspace(ds):
    """Different feature spaces"""
    rc.mincount = 1
    rc.min_infogain = 0
    rc.positives_only = True
    rc.scoremethod = "scores_laplace_split"
    s = base / "FSPACE" / "lsplit_df1_ig0_pos" / ds / ("sd%d"%rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for fs in spaces:
        fspace, rc.class_mask = spaces[fs]
        base_valid(s/fs, ds, fspace, tab)


def compare_all(*dslist):
    """Perform all three aspect-comparisons (prior, featurespace and Document
    Frequency), on each of the given data sets."""
    for ds in dslist:
        compare_prior(ds)
        compare_featselection(ds)
        compare_featspace(ds)


def compare_wordextract(ds):
    """Different ways to extract word features"""
    rc.mincount = 1
    rc.min_infogain = 0
    rc.positives_only = False
    rc.scoremethod = "scores_laplace_split"
    rc.class_mask = []
    s = base / "WORD" / "lsplit_df1_ig0_all" / ds / ("sd%d"%rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for fs in ["word", "word_folded", "word_nodash", "word_num"]:
        base_valid(s/fs, ds, "feats_"+fs, tab)
    base_valid(s/"iedbword", ds, "feats_iedb_word", tab)


def bmc(*datasets):
    """Do the cross-validation on the sample topics from the BMC manuscript."""
    fs = "meshqi"
    rc.scoremethod = "scores_bgfreq"
    rc.mincount = 1
    rc.min_infogain = 0
    rc.positives_only = False
    fspace, cmask = spaces[fs]
    rc.class_mask = cmask
    s = base / "bmc" / "validation"
    for ds in datasets:
        base_valid(s/ds, ds, fspace)


def gdsmall():
    """Basic cross validation sanity test"""
    ds = "gdsmall"
    fs = "meshqi"
    rc.scoremethod = "scores_bgfreq"
    rc.mincount = 1
    rc.min_infogain = 0
    rc.positives_only = False
    fspace, rc.class_mask = spaces[fs]
    s = base / ds / "validation"
    base_valid(s, ds, fspace)


if __name__ == "__main__":
    iofuncs.start_logger(logfile=False)
    if sys.argv[1] in locals():
        adata = ArticleData.Defaults()
    # Call the named function with provided arguments
    locals()[sys.argv[1]](*sys.argv[2:])


'''
def trec():
    """Performs split-sample validation on the 2005 TREC Genomics track
    categorisation subtasks (Allele, Expression, GO, Tumor)."""
    groupdir = rc.root / "results" / "TREC"
    if not groupdir.exists(): groupdir.mkdir()
    fdata = FeatureData.Defaults("feats_word_mesh_all_filt", nx.uint32)
    rc.mincount = 4
    rc.scoremethod = "laplace_split"
    ntest = rc.corpora / "TREC" / "NEG_test.txt"
    ntrain = rc.corpora / "TREC" / "NEG_train.txt"
    for ds, Ur in [("A",17.0), ("E",64.0), ("G",11.0), ("T",231.0)]:
        rc.utility_r = Ur
        dataset = "TREC_%s" % ds
        ptrain = rc.corpora / "TREC" / (ds + "train.txt")
        ptest = rc.corpora / "TREC" / (ds + "test.txt")
        op = SplitValidation(groupdir/dataset, dataset, adata, fdata)
        op.validation(ptrain, ntrain, ptest, ntest)
    fdata.close()
'''
