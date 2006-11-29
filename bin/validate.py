#!env python

"""Calculate performance statistics

CGI script may provide batchid, number of negatives, number of folds,
and pseudocount as parameters.

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
import article
import genedrug
import medline
import validation

def do_validation():
    # Read training data
    statfile = article.StatusFile(c.statfile, c.dataset, c.nfolds)
    featmap = article.FeatureMapping(c.featuremap)
    featdb = medline.FeatureDatabase(c.featuredb, 'r', dbname=c.featureset)
    positives = set(article.readPMIDFile(c.posfile, featdb))
    negatives = set(article.readPMIDFile(c.negfile, featdb))
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
        pos_arts = article.getArticles(c.articledb, c.posfile)
        neg_arts = article.getArticles(c.articledb, c.negfile)
        gdfilter = validation.getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
        #for article in chain(pos_arts, neg_arts):
        for art in pos_arts:
            gdresult = gdfilter(art)
            art.genedrug = gdresult
            if len(gdresult) > 0:
                genedrug_articles.add(art.pmid)
        from dbexport import countGeneDrug
        gdcount = countGeneDrug(pos_arts)
        outf = file("/tmp/output.csv", "w")
        outf.write("PMID,GENE,DRUG\n")
        for (gene,drug),pmids in gdcount.iteritems():
            for pmid in pmids:
                outf.write("%s,%s,%s\n" % (pmid, gene, drug))
        outf.close()
    val = validation.Validator(
        featmap,
        featdb,
        positives,
        negatives,
        c.nfolds,
        c.pseudocount,
        c.dodaniel,
        genedrug_articles,
        c.exclude_feats,
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
        results = val.validate(statfile)
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
        article.chooseRandomLines(c.articlelist, c.negfile, numnegs)
    do_validation()
