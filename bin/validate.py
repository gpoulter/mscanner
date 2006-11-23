#!env python

"""Calculate performance statistics

@author: Graham Poulter
                                   

Perform cross-validation on PMIDs listed in
configuration.validation.(posfile, negfile) outputting a performance
report.

CGI script may provide batchid, number of negatives, nfolds, recall,
and pseudocount as parameters.

"""

import cPickle
from itertools import chain
import logging as log
import os
from path import path
import sys
import time

import configuration as c
from article import createStatusFile, FeatureMapping, getArticles, chooseRandomLines, readPMIDFile
from genedrug import getGeneDrugFilter
from medline import FeatureDatabase
from validation import Validator

def do_validation():
    # Read training data
    meshdb = FeatureMapping(c.meshdb)
    featdb = FeatureDatabase(c.featuredb, 'r', dbname="meshterms")
    positives = set(readPMIDFile(c.posfile, featdb))
    negatives = set(readPMIDFile(c.negfile, featdb))
    # Remove positives found in the negative set
    log.debug("Removing negatives found in positives") 
    for x in positives:
        if x in negatives:
            negatives.remove(x)
    # Get which document ids have gene-drug assocations
    genedrug_articles = None
    if c.dogenedrug:
        log.debug("Getting gene-drug associations") 
        genedrug_articles = set()
        pos_arts = getArticles(c.articledb, c.posfile)
        neg_arts = getArticles(c.articledb, c.negfile)
        gdfilter = getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
        #for article in chain(pos_arts, neg_arts):
        for article in pos_arts:
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
        c.nfolds,
        c.pseudocount,
        c.dodaniel,
        genedrug_articles
        )
    # Load pickled results
    pickle = c.valid_report / "results.pickle"
    if pickle.isfile():
        log.info("Using cached results from %s", pickle)
        results = cPickle.load(file(pickle, "rb"))
    else:
        # Create directory if necessary
        if not c.valid_report.exists():
            c.valid_report.makedirs()
        # Recalculate results
        log.info("Recalculating results, to store in %s", pickle)
        results = val.validate(c.statfile)
        cPickle.dump(results, file(pickle,"wb"), protocol=2)
    # Output performance statistics
    log.debug("Writing performance statistics")
    val.report(results[0], results[1], c.valid_report, c.stylesheet)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        c.dataset = sys.argv[1]
        numnegs = int(sys.argv[2])
        c.nfolds = int(sys.argv[3])
        c.pseudocount = float(sys.argv[4])
        c.valid_report = c.weboutput / c.dataset
        c.posfile = c.valid_report / "positives.txt"
        c.negfile = c.valid_report / "negatives.txt"
        chooseRandomLines(c.articlelist, c.negfile, numnegs)
    try:
        createStatusFile(c.statfile, c.dataset, c.nfolds)
        do_validation()
    finally:
        c.statfile.remove()
