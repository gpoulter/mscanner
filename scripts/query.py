#!/usr/bin/env python

"""Performs queries using datasets from the MScanner paper

Usage::

    python query.py <dataset> [ <dataset> ... ]

Example::

    python query.py pg07 radiology aids
"""


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


dataset_map = {
    "aidsbio"     : "Paper/aidsbio.txt",
    "pg07"        : "Paper/pg07_full.txt",
    "radiology"   : "Paper/radiology.txt",
    "gdsmall"     : "Test/gdsmall.txt",
    "invalid"     : "Test/invalid.txt",
}
"""Mapping from dataset code to input file"""


def do_query(*datasets):
    adata = ArticleData.Defaults()
    #fdata = FeatureData.Defaults_MeSH()
    fdata = FeatureData.Defaults_All()
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        QM = QueryManager(rc.working / "query" / dataset, dataset, limit=500,
                          adata=adata, fdata=fdata, threshold=0, prior=None)
        QM.query(rc.corpora / dataset_map[dataset])
        QM.write_report()
    adata.close()
    fdata.close()

if __name__ == "__main__":
    iofuncs.start_logger()
    do_query(*sys.argv[1:])
    

