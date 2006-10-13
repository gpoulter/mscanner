#!/usr/bin/env python

"""Calculate performance statistics

@author: Graham Poulter
                                   

Perform cross-validation on PMIDs listed in
configuration.validation.(posfile, negfile) outputting a performance
report.

CGI script may provide batchid, number of negatives, nfolds, recall,
and pseudocount as parameters.

"""

import cPickle
import logging as log
import os
import time
import sys
from path import path

from configuration import base, genedrug as gd, medline as m, validation as v
from article import FeatureMapping, getArticles
from genedrug import getGeneDrugFilter
from medline import FeatureDatabase
from validation import Validator

def random_subset(ifname,ofname,N):
    """Choose N random lines from ifname and print them to ofname"""
    from random import randint
    lines = file(ifname, "r").readlines()
    size = len(lines)
    outfile = file(ofname, "w")
    if N > size:
        raise ValueError("N > length of file")
    i = 0
    while i < N:
        i += 1
        r = randint(0,size-1)
        outfile.write(lines[r])
        lines[r] = lines[size-1]
        size = size - 1

def do_validation():
    # Read training data
    meshdb = FeatureMapping(m.meshdb)
    featdb = FeatureDatabase(m.featuredb, 'r')
    positives = set( int(line.split()[0]) for line in file(v.posfile,"r") if line.strip()!="" )
    negatives = set( int(line.split()[0]) for line in file(v.negfile,"r") if line.strip()!="" )
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
        pos_arts = getArticles(v.posfile)
        neg_arts = getArticles(v.negfile)
        gdfilter = getGeneDrugFilter(gd.genedrug, gd.drugtable, gd.gapscore)
        for article in chain(pos_arts, neg_arts):
            gd = gdfilter(article)
            if len(gd) > 0:
                genedrug_articles.add(article.pmid)
    # Initialist validator
    val = Validator(
        meshdb,
        featdb,
        positives,
        negatives,
        v.recall,
        v.nfolds,
        v.pseudocount,
        v.daniel,
        genedrug_articles
        )
    # Load pickled results
    pickle = v.prefix+"results.pickle"
    if pickle.isfile():
        log.info("Using cached results from %s", v.prefix.dirname().relpathto(pickle))
        results = cPickle.load(file(pickle, "rb"))
    else:
        # Recalculate results
        log.info("Recalculating results, to store in %s", v.prefix.dirname().relpathto(pickle))
        results = val.validate()
        cPickle.dump(results, file(pickle,"wb"), protocol=2)
    # Output performance statistics
    log.debug("Writing performance statistics")
    val.report(results, v.prefix, v.stylesheet)

def cgi_invocation():
    try:
        batchid = sys.argv[1]
        v.numnegs = int(sys.argv[2])
        v.nfolds = int(sys.argv[3])
        v.recall = float(sys.argv[4])
        v.pseudocount = float(sys.argv[5])
        v.prefix = (base.weboutput / batchid) + "/"
        v.posfile = v.prefix / "positives.txt"
        v.negfile = v.prefix / "negatives.txt"
        random_subset(v.allpmids, v.negfile, v.numnegs)
    except Exception:
        print __doc__
        sys.exit(1)
    # Lock with PID, batch id and start time
    pidfile = path("/var/run/medscan.pid")
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
