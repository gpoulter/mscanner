#!/usr/bin/env python

"""Performs retrieval test analysis, where a subset of the input is used to
query, and the results are compared against the remainder of the input."""

from __future__ import division

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import itertools
import logging
import numpy as nx
from path import path
import pprint as pp
import sys

from mscanner.configuration import rc
from mscanner.core import iofuncs
from mscanner.core.QueryManager import QueryManager
from mscanner.medline.FeatureData import FeatureData
from mscanner.medline.FeatureStream import DateAsInteger
from mscanner.medline import Shelf


def truepos_vs_rank(results, test):
    """Evaluate true positives as function.

    @param results: Result PubMed IDs in decreasing order of score
    
    @param test: Set of known relevant articles (to be counted as true positives)

    @return: Vector of true positives at each rank (false
    positives = rank - true positives)."""
    assert isinstance(test, set)
    TP_total = nx.zeros(len(results), nx.int32)
    for idx, pmid in enumerate(results):
        TP_total[idx] = TP_total[idx-1] if idx > 0 else 0
        if pmid in test:
            TP_total[idx] += 1
    return TP_total


def plot_rank_performance(fname, textfname, tp_vs_rank, 
                          total_relevant, total_irrelevant):
    """Plot recall and precision versus rank.

    @param fname: File in which to draw the graph
    
    @param textfname: File in which to write the text version of the graph

    @param tp_vs_rank: Vector where the values are the number of true positives
    as a function of (where rank = index+1). 
    
    @param total_relevant: Total number of relevant articles to be found.
    
    @param total_irrelevant: The number of irrelevant articles retrieved by
    the query that was manually filtered to yield the total_relevant. We use
    this to calculate the precision of the PubMed query against which we are
    comparing. 
    """
    logging.debug("Plotting Retrieval curve to %s", fname.basename())
    from Gnuplot import Gnuplot, Data
    g = Gnuplot()
    g.title("Query performance against positive data")
    g.ylabel("Recall & Precision")
    g.xlabel("Rank")
    g("set terminal png")
    g("set output '%s'" % fname)
    ranks = nx.arange(1,len(tp_vs_rank)+1)
    recall = tp_vs_rank / total_relevant
    precision = tp_vs_rank / ranks
    total_query = total_relevant + total_irrelevant
    q_prec = total_relevant/total_query
    q_recall = (q_prec * ranks) / total_relevant
    iofuncs.write_lines(textfname, ["%d, %d, %.3f, %.3f" % x for x in 
        itertools.izip(itertools.count(), tp_vs_rank, recall, precision)],
        "Rank, TP, Recall, Precision")
    logging.info("There are %d relevant and %d irrelevant in Query results", 
                 total_relevant, total_irrelevant)
    logging.info("At end, Query recall is 1.0 and precision is %f", q_prec)
    logging.info("At end, MScanner recall is %f and precision is %f",
                 recall[-1], precision[-1])
    g.plot(Data(ranks, recall, title="Recall", with="lines"),
           Data(ranks, precision, title="Precision", with="lines"),
           Data(ranks, q_recall, title="Query Recall", with="lines"),
           Data([0,total_query], [q_prec,q_prec], title="Query Precision", with="lines"))


def compare_iedb_query():
    """We compare MScanner retrieval to the complex PubMed query
    used to create the gold standard for the IEDB classifier."""
    dataset = "IEDBQuery"
    rc.mincount = 0
    rc.min_infogain = 2e-5
    rc.scoremethod = "scores_laplace_split"
    rc.type_mask = []
    outdir = rc.root/ "results" / "iedb-final-retrieval"
    if not outdir.exists(): outdir.makedirs()
    mindate, maxdate = 20040101, 20041231
    t_mindate, t_maxdate = 20040101, 20041231
    IEDB = rc.corpora / "IEDB" / "Query Data"
    trainpos_file = IEDB / "iedb-pos-pre2004.txt"
    testpos_file = IEDB / "iedb-pos-2004.txt"
    testneg_file = IEDB / "iedb-neg-2004.txt"
    N_testnegs = len(testneg_file.lines())
    testpos = set(iofuncs.read_pmids(testpos_file))
    limit = len(testpos) + N_testnegs
    fdata = FeatureData.Defaults("feats_wmqia")
    artdb = Shelf.open(rc.articles_home/rc.articledb, 'r')
    QM = QueryManager(outdir, dataset, limit, artdb=artdb, fdata=fdata,
                      mindate=mindate, maxdate=maxdate,
                      t_mindate=t_mindate, t_maxdate=t_maxdate)
    QM.query(trainpos_file)
    mscanner_results = [p for s,p in QM.results]
    truepos = truepos_vs_rank(mscanner_results, testpos)
    plot_rank_performance(outdir/"pr_rank.png", outdir/"pr_rank.txt",
                          truepos, len(testpos), N_testnegs)
    QM.write_report(maxreport=1000)


if __name__ == "__main__":
    iofuncs.start_logger()
    eval(sys.argv[1])
