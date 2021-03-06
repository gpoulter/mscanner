#!/usr/bin/env python

"""Performs queries using datasets from the MScanner paper

Example::
    python query.py query pg07 radiology aids
"""

import numpy as nx
import sys

from mscanner.configuration import rc
from mscanner.core import iofuncs
from mscanner.core.QueryManager import QueryManager
from mscanner.medline.FeatureData import FeatureData
from mscanner.medline import Shelf


                                     
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

dataset_map = {
    "aidsbio"     : "Paper/aidsbio.txt",
    "pg07"        : "Paper/pg07_full.txt",
    "radiology"   : "Paper/radiology.txt",
    "gdsmall"     : "Test/gdsmall.txt",
    "invalid"     : "Test/invalid.txt",
    "iedb04"      : "IEDB/iedb-pos-pre2004.txt",
}


def base_query(outdir, dataset, featurespace, limit=500, df=None, ig=None):
    """Configurable query"""
    # Make output directory and calculate title
    if not outdir.parent.exists(): 
        outdir.parent.makedirs()
    title = "_".join(base.relpathto(outdir).splitall()[1:])
    # Set document frequency cutoff
    if df is not None:
        rc.mincount = df
    # Set information gain cutoff
    if ig is not None:
        rc.min_infogain = ig
    fdata = FeatureData.Defaults(featurespace)
    QM = QueryManager(outdir, title, limit, artdb, fdata, threshold=0)
    QM.query(rc.corpora / dataset_map[dataset])
    QM.write_report()
    fdata.close()
    

def final(*datasets):
    """Do queries on the sample topics using final classifer."""
    groupdir = base / "final-query"
    if not base.exists(): base.mkdir()
    fdata = FeatureData.Defaults("feats_wmqia")
    rc.mincount = 0
    rc.min_infogain = 2e-5
    for dataset in datasets:
        QM = QueryManager(groupdir / dataset, dataset, limit=1000,
                          artdb=artdb, fdata=fdata, threshold=0, prior=None)
        QM.query(rc.corpora / dataset_map[dataset])
        QM.write_report()
    fdata.close()


def old(*datasets):
    """Do queries on the sample topics from the BMC manuscript."""
    groupdir = base / "bmc" / "query_mqi"
    if not base.exists(): base.mkdir()
    fdata = FeatureData.Defaults("feats_mesh_qual_issn")
    for dataset in datasets:
        QM = QueryManager(groupdir / dataset, dataset, limit=1000,
                          artdb=artdb, fdata=fdata, threshold=0, prior=None)
        QM.query(rc.corpora / dataset_map[dataset])
        QM.write_report()
    fdata.close()


def yael():
    dataset = "query_wmqia_ig1e-5"
    groupdir = base / "yael"
    if not base.exists(): base.mkdir()
    fdata = FeatureData.Defaults("feats_wmqia")
    rc.min_infogain = 1e-5
    QM = QueryManager(groupdir / dataset, dataset, limit=2000,
                      artdb=artdb, fdata=fdata, threshold=0, prior=None)
    QM.query(rc.corpora / "Yael" / "PharmGKB-2008.03.18.txt")
    QM.write_report()
    fdata.close()

    
    
def gdsmall():
    ds = "gdsmall"
    s = base / ds / "query"
    rc.scoremethod = "scores_bgfreq"
    base_query(s/"mqi_bgfreq_df1", ds, "feats_mesh_qual_issn", df=1)
    #rc.scoremethod = "scores_laplace_split"
    #base_query(s/"wmqia_bgfreq_df4", ds, "feats_wmqia", df=4)


if __name__ == "__main__":
    iofuncs.start_logger()
    if sys.argv[1] in locals():
        artdb = Shelf.open(rc.articles_home/rc.articledb, 'r')
    # Call the named function with provided arguments
    locals()[sys.argv[1]](*sys.argv[2:])
    

