#!/usr/bin/env python

"""Query the MEDLINE database with a set of positive articles

CGI script may call with batchid, pseudocount, result limit and score
threshold as arguments.

                                   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
"""

import cPickle
from itertools import chain
import logging as log
import numpy as nx
import os
from path import path
import sys
import time

from mscanner import configuration
from mscanner import dbshelve
from mscanner.dbexport import exportDefault
from mscanner.featuredb import FeatureDatabase, FeatureStream
from mscanner.featuremap import FeatureMapping
from mscanner.genedrug import getGeneDrugFilter
from mscanner.scoring import FeatureScoreInfo, filterDocuments, writeReport
from mscanner.utils import countFeatures, runMailer, StatusFile, readPMIDs, writePMIDScores

def do_query(configvars):
    # Perform query
    c = configvars
    log.debug("Peforming query for dataset %s", c.dataset)
    statfile = StatusFile(c.statfile, dataset=c.dataset, start=c.timestamp)
    try:
        # Load article information
        featdb = FeatureDatabase(c.featuredb, 'r')
        featmap = FeatureMapping(c.featuremap)
        artdb = dbshelve.open(c.articledb, 'r')
        input_pmids = set(readPMIDs(c.reportdir / c.posname, include=featdb))
        statfile.total = len(featdb)-len(input_pmids)

        # Calculate feature score information
        pos_counts = countFeatures(len(featmap), featdb, input_pmids)
        neg_counts = nx.array(featmap.counts, nx.int32) - pos_counts
        feature_info = FeatureScoreInfo(
            pos_counts,  neg_counts,
            len(input_pmids), len(featdb) - len(input_pmids),
            c.pseudocount, featmap, c.exclude_types, c.dodaniel)

        # Load saved results
        if (c.reportdir / c.index_html).isfile():
            log.info("Loading saved results")
            inputs = readPMIDs(c.reportdir / c.posname, withscores=True)
            results = readPMIDs(c.reportdir / c.result_scores, withscores=True)
        # Recalculate results
        else:
            log.info("Recalculating results")
            # Calculate and write result scores
            featstream = FeatureStream(file(c.featurestream, "rb"))
            results = filterDocuments(featstream, feature_info.scores, input_pmids, c.limit, c.threshold, statfile)
            writePMIDScores(c.reportdir / c.result_scores, results)
            # Calculate and write input scores
            inputs = [ (pmid,nx.sum(feature_info.scores[featdb[pmid]])) for pmid in input_pmids ]
            writePMIDScores(c.reportdir / c.posname, inputs)

        # Write result report
        log.debug("Writing report")
        writeReport(inputs, results, feature_info, c, artdb)

        # Export to database
        if c.outputdb is not None:
            log.debug("Getting gene-drug associations on results")
            gdfilter = getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
            gdarticles = []
            for pmid, score in chain(results, inputs):
                a = artdb[str(pmid)]
                a.genedrug = gdfilter(a)
                if len(a.genedrug) > 0:
                    gdarticles.append(a)
            log.debug("Exporting database")
            exportDefault(c.outputdb, gdarticles)
            
    finally:
        log.debug("Cleaning up")
        featdb.close()
        artdb.close()
        del statfile
        runMailer(c.smtp_server, c.emails_path)
        
def choose_query(configvars, dataset):
    c = configvars
    c.dataset = dataset
    c.reportdir = c.query_output / c.dataset
    if dataset == "pg04":
        pos = "pharmgkb-2004.txt"
    elif dataset == "pg07":
        pos = "pharmgkb-070205.txt"
    elif dataset == "aids":
        pos = "aids-bioethics-Oct06.txt"
    elif dataset == "radiology":
        pos = "daniel-radiology.txt"
    elif dataset == "mscanner":
        pos = "mscanner-bibliography.txt"
    elif dataset == "gdsmall":
        pos = "genedrug-small.txt"
    else:
        raise ValueError("Invalid query dataset " + dataset)
    if not isinstance(pos, path):
        pos = c.corpora / pos
    if not c.reportdir.isdir():
        c.reportdir.mkdir()
        pos.copy(c.reportdir / c.posname)
    
if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Please give dataset code or full parameters")
    elif len(sys.argv) == 2:
        choose_query(configuration, sys.argv[1])
    elif len(sys.argv) > 2:
        c = configuration
        c.dataset = sys.argv[1]
        c.pseudocount = float(sys.argv[2])
        c.limit = int(sys.argv[3])
        c.threshold = float(sys.argv[4])
        c.reportdir = c.query_output / c.dataset
        path(sys.argv[5]).copy(c.reportdir / c.posname)
    do_query(configuration)
