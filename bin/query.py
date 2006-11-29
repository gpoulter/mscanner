#!env python

"""Query the MEDLINE database with a set of positive articles

CGI script may call with batchid, pseudocount, result limit and score
threshold as arguments.

@author: Graham Poulter
                                   

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
import genedrug
import medline
import scoring

def do_query():
    # Perform query
    log.debug("Opening databases")
    log.debug("Peforming query for dataset %s", c.dataset)
    statfile = article.StatusFile(c.statfile, c.dataset)
    featmap = article.FeatureMapping(c.featuremap)
    featdb = medline.FeatureDatabase(c.featuredb, 'r', dbname=c.featureset)
    statfile.update(0, len(featdb))
    artdb = dbshelve.open(c.articledb, 'r')
    posids = set(article.readPMIDFile(c.posfile, featdb))
    pos_counts = article.TermCounts(items=(featdb[d] for d in posids))
    bg_counts = article.TermCounts(c.termcounts)
    neg_counts = bg_counts.subtract(pos_counts)
    term_scores = scoring.getTermScores(pos_counts, neg_counts, c.pseudocount,
                                        featmap=featmap, exclude=c.exclude_feats)
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
        results = scoring.filterDocuments(
            ((k,v) for k,v in featdb.iteritems() if k not in posids),
            term_scores, c.limit, c.threshold, statfile)
        cPickle.dump(results, file(pickle,"wb"), protocol=2)
    # Write result report
    log.debug("Writing report")
    scoring.writeReport(
        scores = results,
        featmap = featmap,
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
        gdfilter = genedrug.getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
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
    do_query()
