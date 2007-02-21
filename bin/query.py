#!env python

"""Query the MEDLINE database with a set of positive articles

CGI script may call with batchid, pseudocount, result limit and score
threshold as arguments.

@author: Graham Poulter
                                   

"""

import cPickle
from itertools import chain
import logging as log
import numpy as nx
import os
from path import path
import sys
import time

import configuration as c
import dbshelve
import dbexport
import article
import genedrug
import medline
import scoring

def do_query():
    # Perform query
    log.debug("Peforming query for dataset %s", c.dataset)
    statfile = article.StatusFile(c.statfile, c.dataset)
    try:
        # Load article information
        featdb = medline.FeatureDatabase(c.featuredb, 'r')
        featmap = article.FeatureMapping(c.featuremap)
        artdb = dbshelve.open(c.articledb, 'r')
        input_pmids = set(article.readPMIDs(c.reportdir/c.posfile, include=featdb))
        statfile.total = len(featdb)-len(input_pmids)

        # Calculate feature score information
        pos_counts = article.countFeatures(len(featmap), featdb, input_pmids)
        neg_counts = nx.array(featmap.counts, nx.int32) - pos_counts
        feature_info = scoring.FeatureScoreInfo(
            pos_counts,  neg_counts,
            len(input_pmids), len(featdb) - len(input_pmids),
            c.pseudocount, featmap, c.exclude_types, c.dodaniel)

        # Load saved results
        if (c.reportdir/c.index_file).isfile():
            log.info("Loading saved results")
            inputs = article.readPMIDs(c.reportdir/c.posfile, withscores=True)
            results = article.readPMIDs(c.reportdir/c.query_results_name, withscores=True)
        # Recalculate results
        else:
            log.info("Recalculating results")
            # Calculate and write result scores
            queryids = ((k,v) for k,v in featdb.iteritems() if int(k) not in input_pmids)
            results = scoring.filterDocuments(queryids, feature_info.scores, c.limit, c.threshold, statfile)
            article.writePMIDScores(c.reportdir/c.query_results_name, results)
            # Calculate and write input scores
            inputs = [ (pmid,nx.sum(feature_info.scores[featdb[pmid]])) for pmid in input_pmids ]
            article.writePMIDScores(c.reportdir/c.posfile, inputs)

        # Write result report
        log.debug("Writing report")
        scoring.writeReport(inputs, results, feature_info, c, artdb)

        # Export to database
        if c.outputdb is not None:
            log.debug("Getting gene-drug associations on results")
            gdfilter = genedrug.getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
            gdarticles = []
            for pmid, score in chain(results, inputs):
                a = artdb[str(pmid)]
                a.genedrug = gdfilter(a)
                if len(a.genedrug) > 0:
                    gdarticles.append(a)
            log.debug("Exporting database")
            dbexport.exportDefault(c.outputdb, gdarticles)
            
    finally:
        featdb.close()
        artdb.close()
        del statfile
        article.runMailer(c.smtp_server, c.mailer)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Please give dataset code or full parameters")
    elif len(sys.argv) == 2:
        c.choose_query(sys.argv[1])
    elif len(sys.argv) > 2:
        c.configure_query(
            dataset = sys.argv[1],
            pseudocount = float(sys.argv[2]),
            limit = int(sys.argv[3]),
            threshold = float(sys.argv[4]),
            pos = path(sys.argv[5])
            )
    do_query()
