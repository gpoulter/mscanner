#!/usr/bin/env python

"""Calculate performance statistics

CGI script may provide batchid, number of negatives, number of folds,
pseudocount, and alpha as parameters.

                                   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
"""

from itertools import chain, izip
import logging as log
import numpy as nx
from path import path
import sys

from mscanner import configuration
from mscanner.dbexport import writeGeneDrugCountsCSV
from mscanner.featuredb import FeatureDatabase
from mscanner.featuremap import FeatureMapping
from mscanner.genedrug import getGeneDrugFilter
from mscanner.utils import getArticles, readPMIDs, StatusFile, runMailer, writePMIDScores
from mscanner import validation

def do_validation(configvars):
    c = configvars
    statfile = StatusFile(c.statfile, dataset=c.dataset, total=c.nfolds, start=c.timestamp)
    try:
        # Load up info about features and articles
        featmap = FeatureMapping(c.featuremap)
        featdb = FeatureDatabase(c.featuredb, 'r')

        # Load already-calculated scores
        if (c.reportdir / c.index_html).isfile():
            log.info("Using cached results for %s", c.dataset)
            positives, pscores = zip(*readPMIDs(c.reportdir / c.posname, withscores=True))
            negatives, nscores = zip(*readPMIDs(c.reportdir / c.negname, withscores=True))
            positives = nx.array(positives, nx.int32)
            pscores = nx.array(pscores, nx.float32)
            negatives = nx.array(negatives, nx.int32)
            nscores = nx.array(nscores, nx.float32)

        # Recalculate scores
        else:
            log.info("Recalculating results for %s", c.dataset)
            log.info("Reading positives")
            positives = list(readPMIDs(c.reportdir / c.posname, include=featdb))
            log.info("Reading negatives")
            negatives = list(readPMIDs(c.reportdir / c.negname, exclude=set(positives)))
            positives = nx.array(positives, nx.int32)
            negatives = nx.array(negatives, nx.int32)
            log.info("Done reading")
            # Get which document ids have gene-drug assocations
            genedrug_articles = None
            if c.dogenedrug:
                log.debug("Getting gene-drug associations") 
                genedrug_articles = set()
                pos_arts = getArticles(c.articledb, c.reportdir / c.posname)
                neg_arts = getArticles(c.articledb, c.reportdir / c.negname)
                gdfilter = getGeneDrugFilter(c.genedrug, c.drugtable, c.gapscore)
                for art in chain(pos_arts, neg_arts):
                    gdresult = gdfilter(art)
                    art.genedrug = gdresult
                    if len(gdresult) > 0:
                        genedrug_articles.add(art.pmid)
                writeGeneDrugCountsCSV(dbexport.countGeneDrug(pos_arts))
            val = validation.Validator(
                featmap = featmap,
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
            writePMIDScores(c.reportdir / c.posname, izip(positives, pscores))
            writePMIDScores(c.reportdir / c.negname, izip(negatives, nscores))

        # Output performance statistics
        log.debug("Writing performance statistics")
        validation.report(positives, negatives, pscores, nscores, featmap, featdb, c)
        
    finally:
        del statfile
        runMailer(c.smtp_server, c.emails_path)

def choose_validation(configvars, dataset):
    c = configvars
    c.dataset = dataset
    c.reportdir = c.valid_output / c.dataset
    # Primary results
    if dataset == "aids-vs-500k":
        pos = "aids-bioethics-Oct06.txt"
        neg = "medline07-500k.txt"
    elif dataset == "pg07-vs-500k":
        pos = "pharmgkb-070205.txt"
        neg = "medline07-500k.txt"
    elif dataset == "radiology-vs-500k":
        pos = "daniel-radiology.txt"
        neg = "medline07-500k.txt"
    elif dataset == "random10k-vs-500k":
        pos = "random10k-06.txt"
        neg = "medline07-500k.txt"
    elif dataset == "aids-vs-500k-noissn":
        pos = "aids-bioethics-Oct06.txt"
        neg = "medline07-500k.txt"
        c.exclude_types = ["issn"]
    # Comparing with Daniel on PG04
    elif dataset == "pg04-vs-30k":
        pos = "pharmgkb-2004.txt"
        neg = "medline07-30k.txt"
    elif dataset == "pg04-vs-30k-dan":
        pos = "pharmgkb-2004.txt"
        neg = "medline07-30k.txt"
        c.exclude_types = ["issn"]
        c.nfolds = 10
        c.dodaniel = True
    elif dataset == "pg04-vs-500k":
        pos = "pharmgkb-2004.txt"
        neg = "medline07-500k.txt"    
    # Comparing with Daniel on PG07
    elif dataset == "pg07-vs-30k":
        pos = "pharmgkb-070205.txt"
        neg = "medline07-30k.txt"
    elif dataset == "pg07-vs-500k-dan":
        pos = "pharmgkb-070205.txt"
        neg = "medline07-500k.txt"
        c.exclude_types = ["issn"]
        c.nfolds = 10
        c.dodaniel = True
    elif dataset == "pg07-vs-500k-noissn":
        pos = "pharmgkb-070205.txt"
        neg = "medline07-500k.txt"
        c.exclude_types = ["issn"]
    # Other experiments
    elif dataset == "mscanner-vs-500k":
        pos = "mscanner-bibliography.txt"
        neg = "medline07-500k.txt"
    elif dataset == "pg07-vs-med07":
        pos = "pharmgkb-070205.txt"
        neg = articlelist
    elif dataset == "gdsmall-vs-sample":
        pos = "genedrug-small.txt"
        neg = c.articlelist
    elif dataset == "gdsmall-vs-sample-dan":
        pos = "genedrug-small.txt"
        neg = c.articlelist
        c.exclude_types=["issn"]
        c.nfolds = 10
        c.dodaniel = True
    else:
        raise ValueError("Invalid validation dataset " + c.dataset)
    if not isinstance(pos, path):
        pos = c.corpora / pos
    if not isinstance(neg, path):
        neg = c.corpora / neg
    if not c.reportdir.isdir():
        c.reportdir.mkdir()
        pos.copy(c.reportdir / c.posname)
        neg.copy(c.reportdir / c.negname)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Please give dataset code")
    elif len(sys.argv) == 2:
        choose_validation(configuration, sys.argv[1])
    elif len(sys.argv) > 2:
        c = configuration
        c.dataset = sys.argv[1]
        c.numnegs = int(sys.argv[2])
        c.nfolds = int(sys.argv[3])
        c.pseudocount = float(sys.argv[4])
        c.alpha = float(sys.argv[5])
        c.reportdir = c.valid_output / c.dataset
        path(sys.argv[6]).copy(c.reportdir / c.posname)
        (c.reportdir / c.negname).write_lines(random.sample(c.articlelist.lines(), c.numnegs))
    do_validation(configuration)
