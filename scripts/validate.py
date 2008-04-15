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
from mscanner.core.ValidationManager import CrossValidation
from mscanner.medline.FeatureData import FeatureData

                                     
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
base = rc.root / "results" / "comparisons"

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

    "wordfold":     ("feats_word_fold", []),
    "wordnum":      ("feats_word_num", []),
    "wordstrip":    ("feats_word_strip", []),
    
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

# Cache of loaded data sets (path->array of PubMed IDs)
pmids_cache = {}

class ResultsTable:
    """Tabulate runtime parameters and performance statistics results.
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
            self.stream.write("Name, AvPrec, BEP, ROC, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\n")
        
    def close(self):
        self.stream.close()
        
    def add_results(self, vmanager):
        """Add a row of results from L{CrossValidation}."""
        v = vmanager.metric_vectors
        results = [vmanager.dataset, v.AvPrec, v.breakeven, v.W]
        fmt = "%s, %.3f, %.3f, %.4f" + (", %.2f"*11) + "\n"
        self.stream.write(fmt % tuple(results + v.precision_11))
        self.stream.flush()


def get_dataset(dataset, fdata):
    """Given a data set name, return (pos, neg) pair of arrays of PubMed IDs for
    positive and negative articles."""
    if dataset not in dataset_map:
        raise ValueError("Invalid Data Set %s" % dataset)
    # Get pair of incomplete paths
    dspaths = list(dataset_map[dataset])
    results = []
    # Loop over the two file names (pos and neg)
    for dspath in dspaths: 
        if isinstance(dspath, str):
            # Complete the path, and cache loaded PubMed IDs for later
            try:
                pmids = pmids_cache[dspath]
            except KeyError:
                logging.debug("Reading PMIDs for %s", dspath)
                pmids, notfound, exclude = \
                     iofuncs.read_pmids_careful(rc.corpora / dspath, fdata.featuredb)
                pmids_cache[dspath] = pmids
            results.append(pmids)
        else:
            # Path was not a string, let ValidationManager sort it out
            results.append(dspath)
    # List of [pos,neg] to tuple (pos,neg)
    return tuple(results)


def base_valid(outdir, dataset, featurespace, tab=None, df=None, ig=None, skip=False):
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
    if skip and outdir.exists(): return
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
    # Cache loaded feature spaces (saves a few seconds)
    if featurespace not in fdata_cache:
        fdata_cache[featurespace] = FeatureData.Defaults(featurespace)
    fdata = fdata_cache[featurespace]
    # Perform cross validation
    op = CrossValidation(outdir, title, fdata)
    op.validation(*get_dataset(dataset, fdata))
    op.report_validation()
    if tab is not None:
        tab.add_results(op)


def testing():
    """Miscellaneous individual tests"""
    rc.positives_only = False
    rc.scoremethod = "scores_laplace_split"
    s = base / "df2_ig20_radiology_sd124_wmqia_lsplit"
    fspace, rc.type_mask = spaces["wmqia"]
    base_valid(s, "radiology", fspace, df=2, ig=2e-5)


def compare_prior(*dslist):
    """Different priors"""
    logging.info("COMPARING PRIORS")
    rc.positives_only = False
    s = base / ("smoothing_%d" % rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for ds in dslist:
        for fs in ["wmqia"]:
            fspace, rc.type_mask = spaces[fs]
            for method in ["bgfreq", "laplace_split", "laplace"]:
                rc.scoremethod = "scores_" + method
                base_valid(s/ds/fs/"df2ig20"/method, ds, fspace, tab, df=2, ig=2e-5, skip=True)
                base_valid(s/ds/fs/"ig0"/method, ds, fspace, tab, df=0, ig=2e-5, skip=True)


def compare_featselection(*dslist):
    """Different feature selection approaches."""
    logging.info("COMPARING FEATURE SELECTION")
    rc.scoremethod = "scores_laplace_split"
    s = base / ("selection_%d" % rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for ds in dslist:
        for fs in ["wmqia"]:
            fspace, rc.type_mask = spaces[fs]
            # Varying: Relative Information Gain
            for ig in [1e-6, 1e-5, 2e-5, 1e-4, 1e-3]:
                base_valid(s/ds/fs/("ig%d"%(int(ig*1e6))), ds, fspace, tab, df=0, ig=ig, skip=True)
            # Varying: Document Frequency
            for df in [0, 1, 2, 3, 4, 8]:
                base_valid(s/ds/fs/("df%d"%df), ds, fspace, tab, df=df, ig=0, skip=True)
            # Varying: Whether to use relevant-only features
            rc.positives_only = True
            base_valid(s/ds/fs/"posonly", ds, fspace, tab, df=0, ig=0, skip=True)
            rc.positives_only = False
            # Trying with best combination
            df = 2
            ig = 2e-5
            base_valid(s/ds/fs/("df%dig%d"%(df,int(ig*1e6))), ds, fspace, tab, df=df, ig=ig, skip=True)


def compare_featspace(*dslist):
    """Different feature spaces"""
    logging.info("COMPARING FEATURE SPACES")
    rc.mincount = 0
    rc.min_infogain = 0
    rc.positives_only = False
    rc.scoremethod = "scores_laplace_split"
    s = base / ("featspaces_%d" % rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for ds in dslist:
        for fs in ["meshnq","qual","issn","word","author","meshqi","wmqia","wmqiafilt","iedbconcat"]:
            fspace, rc.type_mask = spaces[fs]
            base_valid(s/ds/"ig0"/fs, ds, fspace, tab, df=0, ig=0, skip=True)
            base_valid(s/ds/"ig20"/fs, ds, fspace, tab, df=0, ig=2e-5, skip=True)


def compare_wordextract(*dslist):
    """Compare ways of extracting word features"""
    logging.info("COMPARING WORD EXTRACTION")
    rc.mincount = 0
    rc.min_infogain = 0
    rc.positives_only = False
    rc.scoremethod = "scores_laplace_split"
    rc.type_mask = []
    s = base / ("wordextract_%d" % rc.randseed)
    tab = ResultsTable(s/"results.txt", append=True)
    for ds in dslist:
        for fs in ["word", "wordnum", "wordfold", "wordstrip", "iedbword" ]:
            fspace, rc.type_mask = spaces[fs]
            base_valid(s/ds/fs, ds, fspace, tab, skip=True)


def compare_all(*dslist):
    """Perform all three aspect-comparisons (prior, featurespace and Document
    Frequency), on each of the given data sets."""
    compare_featselection(*dslist)
    compare_prior(*dslist)
    compare_featspace(*dslist)
    compare_wordextract(*dslist)


def bmc(*datasets):
    """Do the cross-validation on the sample topics from the BMC manuscript."""
    fs = "meshqi"
    #rc.scoremethod = "scores_bgfreq"
    rc.scoremethod = "scores_laplace_split"
    rc.mincount = 1
    rc.min_infogain = 1e-4
    rc.positives_only = False
    fspace, rc.type_mask = spaces[fs]
    s = base / "bmc" / "valid"
    tab = ResultsTable(s/"summary.txt", append=True)
    for ds in datasets:
        base_valid(s/ds, ds, fspace, tab)


def gdsmall():
    """Basic cross validation sanity test"""
    ds = "gdsmall"
    fs = "meshqi"
    rc.scoremethod = "scores_bgfreq"
    rc.mincount = 1
    rc.min_infogain = 0
    rc.positives_only = False
    fspace, rc.type_mask = spaces[fs]
    s = base / ds / "validation"
    base_valid(s, ds, fspace)


if __name__ == "__main__":
    iofuncs.start_logger()
    # Call the named function with provided arguments
    locals()[sys.argv[1]](*sys.argv[2:])
