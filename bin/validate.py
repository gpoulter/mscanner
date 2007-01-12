#!env python

"""Calculate performance statistics

CGI script may provide batchid, number of negatives, number of folds,
and pseudocount as parameters.

@author: Graham Poulter
                                   

"""

#Builtin
import cPickle
from itertools import chain
import logging as log
import os
import random
import sys
import time
#Local
from path import path
import numpy
#MScanner
import configuration as c
import article
import dbexport
import genedrug
import medline
import validation

def do_validation():
    import article
    statfile = article.StatusFile(c.statfile, c.dataset, c.nfolds)
    try:
        featmap = article.FeatureMapping(c.featuremap)
        featdb = medline.FeatureDatabase(c.featuredb, 'r')
        positives = set(article.readPMIDFile(c.posfile, featdb))
        print "READING NEGATIVES"
        negatives = [x for x in article.readPMIDFile(c.negfile, featdb) if x not in positives]
        print "READING NEGATIVES DONE"
        negatives = numpy.array(negatives, dtype=numpy.int32)
        positives = numpy.array(list(positives), dtype=numpy.int32)
        # Get which document ids have gene-drug assocations
        genedrug_articles = None
        if c.dogenedrug:
            log.debug("Getting gene-drug associations") 
            genedrug_articles = set()
            pos_arts = article.getArticles(c.articledb, c.posfile)
            neg_arts = article.getArticles(c.articledb, c.negfile)
            gdfilter = genedrug.getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
            for article in chain(pos_arts, neg_arts):
                gdresult = gdfilter(art)
                art.genedrug = gdresult
                if len(gdresult) > 0:
                    genedrug_articles.add(art.pmid)
            dbexport.writeGeneDrugCountsCSV(dbexport.countGeneDrug(pos_arts))
        val = validation.Validator(
            featmap,
            featdb,
            positives,
            negatives,
            c.nfolds,
            c.pseudocount,
            c.dodaniel,
            genedrug_articles,
            c.exclude_feats,)
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
    finally:
        del statfile

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Please give dataset code")
    elif len(sys.argv) == 2:
        c.choose_validation(sys.argv[1])
    elif len(sys.argv) > 2:
        c.dataset = sys.argv[1]
        numnegs = int(sys.argv[2])
        c.nfolds = int(sys.argv[3])
        c.pseudocount = float(sys.argv[4])
        c.valid_report = c.weboutput / c.dataset
        c.posfile = c.valid_report / "positives.txt"
        c.negfile = c.valid_report / "negatives.txt"
        c.negfile.write_lines(random.sample(c.articlelist.lines(), numnegs))
    do_validation()
