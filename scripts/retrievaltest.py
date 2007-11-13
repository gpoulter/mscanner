#!/usr/bin/env python

"""Performs retrieval test analysis, where a subset of the input is used to
query, and the results are compared against the remainder of the input.

Usage::

    python retrievaltest.py <dataset> <filename>
    
For example::

    python retrievaltest.py pg07 'pharmgkb-070205.txt'
"""

from __future__ import division

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import logging as log
import numpy as nx
from path import path
import pprint as pp
import random
import sys

from mscanner.configuration import rc, start_logger
from mscanner.core import iofuncs
from mscanner.core.QueryManager import QueryManager
from mscanner.medline.Databases import Databases
from mscanner.medline.FeatureStream import Date2Integer


def result_relevance(results, test):
    """Test of query results against a gold standard.

    @param results: List of result PubMed IDs in decreasing order of score
    
    @param test: Set of known relevant articles that are to be counted
    as true positive in L{results}

    @return: Vector storing number of true positives at each rank. False
    positives can be calculated as rank minus true positives."""
    assert isinstance(test, set)
    TP_total = nx.zeros(len(results), nx.int32)
    for idx, pmid in enumerate(results):
        TP_total[idx] = TP_total[idx-1] if idx > 0 else 0
        if pmid in test:
            TP_total[idx] += 1
    return TP_total


def plot_rank_performance(fname, tp_vs_rank, total_relevant, compare_irrelevant):
    """Plot recall and precision versus rank.
    @param fname: File to write graph to

    @param tp_vs_rank: Vector where the values are the number of true positives
    at rank equal to current index plus one. 
    
    @param total_relevant: Total number of relevant articles to be found.
    
    @param compare_irrelevant: The number of irrelevant articles retrieved by
    the query that was manually filtered to yield the L{total_relevant}. We use
    this to calculate the precision of the PubMed query against which we are
    comparing. """
    log.debug("Plotting Retrieval curve to %s", fname)
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
    #print "Ranks: ", ranks
    #print "Recall: ", recall
    #print "Precision: ", precision
    total_query = total_relevant + compare_irrelevant
    g.plot(Data(ranks, recall, title="Recall", with="lines"),
           Data(ranks, precision, title="Precision", with="lines"),
           Data([total_query], [1.0], title="QRecall", with="points pt 3 ps 1"),
           Data([total_query], [total_relevant/total_query], title="QPrecision", with="points pt 5 ps 1"))



def compare_pg07_and_pgxquery():
    """Carry out retrieval comparison for a query to MScanner
    with PG07 versus a query to PubMed with 'pharmacogenetics'
    
    The relevant examples are split randomly into train/test.  A query is done
    with train, and we evaluate the number of members of test that show
    up in the query results as a function of rank"""
    dataset = "pg07_retrieve"
    pos_file = rc.corpora / "Paper" / "pharmgkb_2007.02.05.txt"
    outdir = rc.working / "query" / dataset
    test_proportion = 0.2
    # Load and split PubMed IDs into train/test
    QM = QueryManager(outdir, dataset, limit=1000)
    pmids_list = list(iofuncs.read_pmids(pos_file))
    random.shuffle(pmids_list)
    subset_size = int(test_proportion * len(pmids_list))
    train_pmids = pmids_list[:subset_size]
    test_pmids = pmids_list[subset_size:]
    iofuncs.write_lines(outdir/"test_pmids.txt", test_pmids)
    QM.query(train_pmids)
    QM.env.close()
    # Evaluate true positives as a function of rank
    tp_vs_rank = result_relevance(
        [p for s,p in QM.results], test_pmids)
    iofuncs.write_lines(outdir/"tp_vs_rank.txt", enumerate(tp_vs_rank))



def split_by_date(artdb, pmids, splitdate, before_file, after_file):
    """Split a list of PubMed IDs by record date and save them to file. It
    reloads the results instead if the files already exist.
    
    @param artdb: Mapping from PubMed ID to Article object.
    @param pmids: Sequence of PubMed IDs.
    @param splitdate: Partitions are "before" and "on-or-after" splitdate
    @return: before, after lists of PubMed IDs
    """
    if before_file.exists() and after_file.exists():
        log.info("Loading PMIDs from %s and %s", 
                 before_file.basename(), after_file.basename())
        before = list(iofuncs.read_pmids(before_file))
        after = list(iofuncs.read_pmids(after_file))
    else:
        before, beforedates = [], []
        after, afterdates = [], []
        for pmid in pmids:
            if str(pmid) in artdb:
                date = Date2Integer(artdb[str(pmid)].date_completed)
                if date < splitdate:
                    before.append(pmid)
                    beforedates.append(date)
                else:
                    after.append(pmid)
                    afterdates.append(date)
            else:
                log.error("Failed to find %d.", pmid)
        iofuncs.write_lines(
            before_file, sorted(zip(before, beforedates), key=lambda x:x[1]))
        iofuncs.write_lines(
            after_file, sorted(zip(after, afterdates), key=lambda x:x[1]))
    return before, after


def compare_wang2007_query(test=False):
    """We compare MScanner retrieval to the complex PubMed query
    used in the Wang2007 paper."""
    dataset = "wang2007query"
    outdir = rc.working / "query" / dataset
    if not outdir.exists(): outdir.makedirs()
    if test:
        splitdate = 19920101
        limit = 10000
        pos_file = rc.corpora / "Test" / "gdsmall.txt"
        neg_file = rc.articlelist
    else:
        splitdate = 20050101
        limit = 10000
        pos_file = rc.corpora / "Wang2007" / "combined_pos.txt"
        neg_file = rc.corpora / "Wang2007" / "combined_neg.txt"
    pos_pmids = list(iofuncs.read_pmids(pos_file))
    neg_pmids = list(iofuncs.read_pmids(neg_file))
    env = Databases()
    old_pos, new_pos = split_by_date(env.artdb, pos_pmids, splitdate,
                                     outdir/"old_pos.txt", outdir/"new_pos.txt")
    old_neg, new_neg = split_by_date(env.artdb, neg_pmids, splitdate,
                                     outdir/"old_neg.txt", outdir/"new_neg.txt")
    QM = QueryManager(outdir, dataset, limit, mindate=splitdate, env=env)
    QM.query(old_pos)
    new_results = [p for s,p in QM.results]
    truepos = result_relevance(new_results, set(new_pos))
    iofuncs.write_lines(outdir/"tp_vs_rank.txt", enumerate(truepos))
    plot_rank_performance(outdir/"perf_vs_rank.png", 
                          truepos, len(new_pos),
                          compare_irrelevant=len(new_neg))
    QM.write_report(maxreport=1000)
    env.close()
    


if __name__ == "__main__":
    start_logger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression to execute"
    else:
        eval(sys.argv[1])
