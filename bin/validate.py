#!env python

"""Calculate performance statistics

@author: Graham Poulter
                                   

Perform cross-validation on PMIDs listed in
configuration.validation.(posfile, negfile) outputting a performance
report.

CGI script may provide batchid, number of negatives, nfolds, recall,
and pseudocount as parameters.

"""

from itertools import chain
import cPickle
import logging as log
import os
import time
import sys
from path import path

from configuration import base, genedrug as gd, medline as m, validation as v
from article import FeatureMapping, getArticles, chooseRandomLines, readPMIDFile
from genedrug import getGeneDrugFilter
from medline import FeatureDatabase
from validation import Validator

def do_validation():
    # Read training data
    meshdb = FeatureMapping(m.meshdb)
    featdb = FeatureDatabase(m.featuredb, 'r')
    positives = set(readPMIDFile(v.posfile, featdb))
    negatives = set(readPMIDFile(v.negfile, featdb))
    # Remove positives found in the negative set
    log.debug("Removing negatives found in positives") 
    for x in positives:
        if x in negatives:
            negatives.remove(x)
    # Get which document ids have gene-drug assocations
    genedrug_articles = None
    if v.genedrug:
        log.debug("Getting gene-drug associations") 
        genedrug_articles = set()
        pos_arts = getArticles(m.articledb, v.posfile)
        neg_arts = getArticles(m.articledb, v.negfile)
        gdfilter = getGeneDrugFilter(gd.genedrug, gd.drugtable, gd.gapscore)
        for article in chain(pos_arts, neg_arts):
            gdresult = gdfilter(article)
            article.genedrug = gdresult
            if len(gdresult) > 0:
                genedrug_articles.add(article.pmid)
        from dbexport import countGeneDrug
        gdcount = countGeneDrug(pos_arts)
        outf = file("/tmp/output.csv", "w")
        outf.write("PMID,GENE,DRUG\n")
        for (gene,drug),pmids in gdcount.iteritems():
            for pmid in pmids:
                outf.write("%s,%s,%s\n" % (pmid, gene, drug))
        outf.close()
    # Initialist validator
    val = Validator(
        meshdb,
        featdb,
        positives,
        negatives,
        v.nfolds,
        v.pseudocount,
        v.daniel,
        genedrug_articles
        )
    # Load pickled results
    pickle = v.prefix / "results.pickle"
    if pickle.isfile():
        log.info("Using cached results from %s", pickle)
        results = cPickle.load(file(pickle, "rb"))
    else:
        # Create directory if necessary
        if not v.prefix.exists():
            v.prefix.makedirs()
        # Recalculate results
        log.info("Recalculating results, to store in %s", pickle)
        results = val.validate()
        cPickle.dump(results, file(pickle,"wb"), protocol=2)
    # Output performance statistics
    log.debug("Writing performance statistics")
    val.report(results[0], results[1], v.prefix, v.stylesheet)

def cgi_invocation():
    batchid = sys.argv[1]
    v.numnegs = int(sys.argv[2])
    v.nfolds = int(sys.argv[3])
    v.pseudocount = float(sys.argv[4])
    v.prefix = base.weboutput / batchid
    v.posfile = v.prefix / "positives.txt"
    v.negfile = v.prefix / "negatives.txt"
    chooseRandomLines(m.articlelist, v.negfile, v.numnegs)
    # Lock with PID, batch id and start time
    pidfile = path("/var/run/mscanner.pid")
    if pidfile.exists():
        raise RuntimeError("MedScanner instance is already running")
    pidfile.write_text("%d\n%s\n%s\n" % (os.getpid(), batchid, str(time.time())))
    try:
        do_validation()
    finally:
        pidfile.remove()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cgi_invocation()
    else:
        do_validation()
