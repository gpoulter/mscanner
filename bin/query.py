#!env python

"""Query the MEDLINE database with a set of positive articles

CGI script may call with batchid, pseudocount, result limit and score
threshold as arguments.

@author: Graham Poulter
                                   

"""

#Builtin
import cPickle
from itertools import chain
import logging as log
import os
import sys
import time
#Local
from path import path
import numpy
#MScanner
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
        featmap = article.FeatureMapping(c.featuremap)
        featdb = medline.FeatureDatabase(c.featuredb, 'r')
        statfile.update(0, len(featdb))
        artdb = dbshelve.open(c.articledb, 'r')
        posids = set(article.readPMIDFile(c.posfile, featdb)) # set of int
        pfreqs = article.countFeatures(len(featmap), featdb, posids)
        nfreqs = featmap.featureCounts() - pfreqs
        pdocs = len(posids)
        ndocs = len(featdb) - len(posids)
        termscores = scoring.getTermScores(
            pfreqs, nfreqs, pdocs, ndocs,
            c.pseudocount, featmap=featmap, exclude=c.exclude_feats)
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
            queryids = ((k,v) for k,v in featdb.iteritems() if int(k) not in posids)
            results = scoring.filterDocuments(queryids,         
                termscores, c.limit, c.threshold, statfile)
            cPickle.dump(results, file(pickle,"wb"), protocol=2)
        # Write result report
        log.debug("Writing report")
        scoring.writeReport(
            scores = results,
            featmap = featmap,
            pdocs = pdocs,
            ndocs = ndocs,
            termscores = termscores,
            prefix = c.query_report,
            stylesheet = c.stylesheet,
            pseudocount = c.pseudocount,
            limit = c.limit,
            threshold = c.threshold,
            dataset = c.dataset,
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
    finally:
        del statfile
        article.runMailer(c.smtp_server, c.mailer)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Please give dataset code")
    elif len(sys.argv) == 2:
        c.choose_query(sys.argv[1])
    elif len(sys.argv) > 2:
        c.dataset = sys.argv[1]
        c.pseudocount = float(sys.argv[2])
        c.limit = int(sys.argv[3])
        c.threshold = float(sys.argv[4])
        c.query_report = (c.weboutput / c.dataset) + "/"
        c.posfile = c.query_report / "positives.txt"
    do_query()
