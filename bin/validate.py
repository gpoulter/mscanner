#!env python

"""Calculate performance statistics

CGI script may provide batchid, number of negatives, number of folds,
 pseudocount, and alpha as parameters.

@author: Graham Poulter
                                   

"""

#Builtin
import cPickle
from itertools import chain, izip
import logging as log
import os
import sys
#Local
from path import path
import numpy
#MScanner
import configuration as c
import article
import dbexport
import genedrug
import medline
import scoring
import validation

def do_validation():
    statfile = article.StatusFile(c.statfile, c.dataset, c.nfolds)
    try:
        # Load up info about features and articles
        featmap = article.FeatureMapping(c.featuremap)
        featdb = medline.FeatureDatabase(c.featuredb, 'r')

        # Load already-calculated scores
        if (c.reportdir/c.index_file).isfile():
            log.info("Using cached results for %s", c.dataset)
            positives, pscores = zip(*article.readPMIDs(c.reportdir/c.posfile, withscores=True))
            negatives, nscores = zip(*article.readPMIDs(c.reportdir/c.negfile, withscores=True))
            positives = numpy.array(positives, numpy.int32)
            pscores = numpy.array(pscores, numpy.float32)
            negatives = numpy.array(negatives, numpy.int32)
            nscores = numpy.array(nscores, numpy.float32)

        # Recalculate scores
        else:
            log.info("Recalculating results for %s", c.dataset)
            log.info("Reading positives")
            positives = list(article.readPMIDs(c.reportdir/c.posfile, include=featdb))
            log.info("Reading negatives")
            negatives = list(article.readPMIDs(c.reportdir/c.negfile, exclude=set(positives)))
            positives = numpy.array(positives, numpy.int32)
            negatives = numpy.array(negatives, numpy.int32)
            log.info("Done reading")
            # Get which document ids have gene-drug assocations
            genedrug_articles = None
            if c.dogenedrug:
                log.debug("Getting gene-drug associations") 
                genedrug_articles = set()
                pos_arts = article.getArticles(c.articledb, c.reportdir/c.posfile)
                neg_arts = article.getArticles(c.articledb, c.reportdir/c.negfile)
                gdfilter = genedrug.getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
                for art in chain(pos_arts, neg_arts):
                    gdresult = gdfilter(art)
                    art.genedrug = gdresult
                    if len(gdresult) > 0:
                        genedrug_articles.add(art.pmid)
                dbexport.writeGeneDrugCountsCSV(dbexport.countGeneDrug(pos_arts))
            val = validation.Validator(
                numfeats = len(featmap),
                featdb = featdb,
                pos = positives,
                neg = negatives,
                nfold = c.nfolds,
                pseudocount = c.pseudocount,
                daniel = c.dodaniel,
                genedrug_articles = genedrug_articles,
                mask = featmap.featureTypeMask(c.exclude_types)
                )
            pscores, nscores = val.validate(statfile)
            article.writePMIDScores(c.reportdir/c.posfile, izip(positives, pscores))
            article.writePMIDScores(c.reportdir/c.negfile, izip(negatives, nscores))

        # Output performance statistics
        log.debug("Writing performance statistics")
        validation.report(positives, negatives, pscores, nscores, featmap, featdb, c)
        
    finally:
        del statfile
        article.runMailer(c.smtp_server, c.mailer)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Please give dataset code")
    elif len(sys.argv) == 2:
        c.choose_validation(sys.argv[1])
    elif len(sys.argv) > 2:
        c.configure_validation(
            dataset = sys.argv[1],
            numnegs = int(sys.argv[2]),
            nfolds = int(sys.argv[3]),
            pseudocount = float(sys.argv[4]),
            alpha = float(sys.argv[5]),
            pos = path(sys.argv[6]),
            )
    do_validation()
