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
import os
import time
import sys
from itertools import chain
import logging as log
from path import path

from configuration import base, genedrug as gd, query as q, medline as m
import dbshelve
import dbexport
from article import FeatureMapping, TermCounts, readPMIDFile
from genedrug import getGeneDrugFilter
from medline import FeatureDatabase
from scoring import getTermScores, filterDocuments, writeReport

def do_query():
    # Perform query
    log.debug("Opening databases")
    meshdb = FeatureMapping(m.meshdb)
    featdb = FeatureDatabase(m.featuredb,'r')
    artdb = dbshelve.open(m.articledb,'r')
    posids = set(readPMIDFile(q.posfile, featdb))
    pos_counts = TermCounts(featdb[d] for d in posids)
    bg_counts = TermCounts.load(m.termcounts)
    neg_counts = bg_counts.subtract(pos_counts)
    term_scores = getTermScores(pos_counts, neg_counts, q.pseudocount)
    # Load pickled results
    pickle = q.prefix / "results.pickle"
    if pickle.isfile():
        log.info("Loading pickled results from %s", pickle)
        results = cPickle.load(file(pickle, "rb"))
    # Recalculate results
    else:
        # Create directory if necessary
        if not q.prefix.exists():
            q.prefix.makedirs()
        log.info("Recalculating results, to store in %s", pickle)
        results = filterDocuments(((k,v) for k,v in featdb.iteritems() if k not in posids),
                                  term_scores, q.limit, q.threshold)
        cPickle.dump(results, file(pickle,"wb"), protocol=2)
    # Write result report
    log.debug("Writing report")
    writeReport(
        scores = results,
        meshdb = meshdb,
        termscores = term_scores,
        pfreq = pos_counts,
        nfreq = neg_counts,
        prefix = q.prefix,
        stylesheet = q.stylesheet,
        pseudocount = q.pseudocount,
        limit = q.limit,
        posfile = q.posfile,
        articles = artdb,
        )
    # Database export
    if q.outputdb is not None:
        log.debug("Getting gene-drug associations on results")
        gdfilter = getGeneDrugFilter(gd.genedrug, gd.drugtable, gd.gapscore)
        gdarticles=[]
        for pmid in chain((a for s,a in results),posids):
            a = artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                gdarticles.append(a)
        log.debug("Exporting database")
        dbexport.exportDefault(q.outputdb, gdarticles)
    featdb.close()
    artdb.close()

def cgi_invocation():
    batchid = sys.argv[1]
    q.pseudocount = float(sys.argv[2])
    q.limit = int(sys.argv[3])
    q.threshold = float(sys.argv[4])
    q.prefix = (base.weboutput / batchid) + "/"
    q.posfile = q.prefix / "positives.txt"
    pidfile = path("/var/run/mscanner.pid")
    if pidfile.exists():
        raise RuntimeError("MedScanner instance is already running")
    # Lock with PID, batch id and start time
    pidfile.write_text("%d\n%s\n%s\n" % (os.getpid(), batchid, str(time.time())))
    try:
        do_query()
    finally:
        pidfile.remove()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cgi_invocation()
    else:
        do_query()
