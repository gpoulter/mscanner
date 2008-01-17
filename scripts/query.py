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
    # Choose feature types
    if featurespace == "feats_mesh_qual_issn":
        ftype = nx.uint16
    else:
        ftype = nx.uint32
    fdata = FeatureData.Defaults(featurespace, ftype)
    QM = QueryManager(outdir, title, limit, adata, fdata, threshold=0)
    QM.query(rc.corpora / dataset_map[dataset])
    QM.write_report()
    fdata.close()
    

def bmc(*datasets):
    """Do queries on the sample topics from the BMC manuscript."""
    groupdir = base / "bmc" / "query"
    if not base.exists(): base.mkdir()
    fdata = FeatureData.Defaults("feats_mesh_qual_issn", nx.uint16)
    for dataset in datasets:
        QM = QueryManager(groupdir / dataset, dataset, limit=1000,
                          adata=adata, fdata=fdata, threshold=0, prior=None)
        QM.query(rc.corpora / dataset_map[dataset])
        QM.write_report()
    fdata.close()
    
    
def gdsmall():
    ds = "gdsmall"
    s = base / ds / "query"
    rc.scoremethod = "scores_bgfreq"
    base_query(s/"MQI_bgfreq_df0", ds, "feats_mesh_qual_issn", df=0)
    base_query(s/"WMQIA_bgfreq_df4", ds, "feats_word_mqi_all", df=4)
    rc.scoremethod = "scores_laplace_split"
    base_query(s/"WMQIF_lsplit_df4", ds, "feats_word_mqi_filt", df=4)


if __name__ == "__main__":
    iofuncs.start_logger(logfile=False)
    if sys.argv[1] in locals():
        adata = ArticleData.Defaults()
    # Call the named function with provided arguments
    locals()[sys.argv[1]](*sys.argv[2:])
    

