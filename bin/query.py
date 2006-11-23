#!env python

"""Query the MEDLINE database with a set of positive articles

@author: Graham Poulter
                                   

Usage: ./query.py

Query the MEDLINE database by training a classifier on
configuration.query.posfile, and output results to configured
locations.

CGI script may call this with batchid, pseudocount and limit as
arguments.

"""

import cPickle
from itertools import chain
import logging as log
import os
from path import path
import sys
import time

import configuration as c
import dbshelve
import dbexport
import article
from genedrug import getGeneDrugFilter
from medline import FeatureDatabase
from scoring import getTermScores, filterDocuments, writeReport

def do_query():
    # Perform query
    log.debug("Opening databases")
    meshdb = article.FeatureMapping(c.meshdb)
    featdb = FeatureDatabase(c.featuredb, 'r', dbname="meshterms")
    article.updateStatusFile(c.statfile, 0, len(featdb))
    artdb = dbshelve.open(c.articledb, 'r')
    posids = set(article.readPMIDFile(c.posfile, featdb))
    pos_counts = article.TermCounts(featdb[d] for d in posids)
    bg_counts = article.TermCounts.load(c.termcounts)
    neg_counts = bg_counts.subtract(pos_counts)
    term_scores = getTermScores(pos_counts, neg_counts, c.pseudocount)
    # Load pickled results
    pickle = c.query_report / "results.pickle"
    if pickle.isfile():
        log.info("Loading pickled results from %s", pickle)
        results = cPickle.load(file(pickle, "rb"))
    # Recalculate results
    else:
        # Create directory if necessary
        if not c.query_report.exists():
            c.query_report.makedirs()
        log.info("Recalculating results, to store in %s", pickle)
        results = filterDocuments(((k,v) for k,v in featdb.iteritems() if k not in posids),
                                  term_scores, c.limit, c.threshold, c.statfile)
        cPickle.dump(results, file(pickle,"wb"), protocol=2)
    # Write result report
    log.debug("Writing report")
    writeReport(
        scores = results,
        meshdb = meshdb,
        termscores = term_scores,
        pfreq = pos_counts,
        nfreq = neg_counts,
        prefix = c.query_report,
        stylesheet = c.stylesheet,
        pseudocount = c.pseudocount,
        limit = c.limit,
        posfile = c.posfile,
        articles = artdb,
        )
    # Database export
    if c.outputdb is not None:
        log.debug("Getting gene-drug associations on results")
        gdfilter = getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
        gdarticles = []
        for pmid in chain((a for s,a in results),posids):
            a = artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                gdarticles.append(a)
        log.debug("Exporting database")
        dbexport.exportDefault(c.outputdb, gdarticles)
    featdb.close()
    artdb.close()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        c.dataset = sys.argv[1]
        c.pseudocount = float(sys.argv[2])
        c.limit = int(sys.argv[3])
        c.threshold = float(sys.argv[4])
        c.query_report = (c.weboutput / c.dataset) + "/"
        c.posfile = c.query_report / "positives.txt"
    try:
        article.createStatusFile(c.statfile, c.dataset)
        do_query()
    finally:
        c.statfile.remove()
