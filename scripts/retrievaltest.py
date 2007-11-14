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

import itertools
import logging
import numpy as nx
from path import path
import pprint as pp
import random
import sys

from mscanner.configuration import rc
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


def plot_rank_performance(fname, textfname, tp_vs_rank, total_relevant, compare_irrelevant):
    """Plot recall and precision versus rank.
    @param fname: File in which to draw the graph
    @param textfname: File in which to write the text version of the graph

    @param tp_vs_rank: Vector where the values are the number of true positives
    at rank equal to current index plus one. 
    
    @param total_relevant: Total number of relevant articles to be found.
    
    @param compare_irrelevant: The number of irrelevant articles retrieved by
    the query that was manually filtered to yield the L{total_relevant}. We use
    this to calculate the precision of the PubMed query against which we are
    comparing. """
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
    total_query = total_relevant + compare_irrelevant
    q_prec = total_relevant/total_query
    q_recall = (q_prec * ranks) / total_relevant
    iofuncs.write_lines(textfname, zip(
        itertools.count(), tp_vs_rank, recall, q_recall, precision),
        "rank\tTP\trecall\tq_recall\tprecision")
    logging.info("There are %d relevant and %d irrelevant in query results", 
                 total_relevant, compare_irrelevant)
    logging.info("At end, Query recall is 1.0 and precision is %f", q_prec)
    logging.info("At end, MScanner recall is %f and precision is %f",
                 recall[0], precision[-1])
    g.plot(Data(ranks, recall, title="Recall", with="lines"),
           Data(ranks, precision, title="Precision", with="lines"),
           Data(ranks, q_recall, title="Query Recall", with="lines"),
           Data([0,total_query], [q_prec,q_prec], title="Query Precision", with="lines"))



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



def split_by_date(artdb, pmids, mindate, maxdate, before_file, after_file, notfound_file):
    """Split a list of PubMed IDs by record date and save them to file. It
    reloads the results instead if the files already exist.
    
    @param artdb: Mapping from PubMed ID to Article object.
    @param pmids: Sequence of PubMed IDs.
    @param mindate: Partitions are "before" and "on-or-after" splitdate
    @param maxdate: Ignore PubMed IDs added after maxdate
    @return: before, after lists of PubMed IDs
    """
    if before_file.exists() and after_file.exists():
        logging.info("Loading PMIDs from %s and %s", 
                 before_file.basename(), after_file.basename())
        before = list(iofuncs.read_pmids(before_file))
        after = list(iofuncs.read_pmids(after_file))
        notfound = list(iofuncs.read_pmids(notfound_file))
    else:
        before, beforedates = [], []
        after, afterdates = [], []
        notfound, ignore = [], []
        for pmid in pmids:
            if str(pmid) in artdb:
                date = Date2Integer(artdb[str(pmid)].date_completed)
                if date < mindate:
                    before.append(pmid)
                    beforedates.append(date)
                elif date <= maxdate:
                    after.append(pmid)
                    afterdates.append(date)
                else:
                    ignore.append(pmid)
            else:
                notfound.append(pmid)
        logging.info("PMIDs: %d before, %d inrange, %d outrange, %d unknown.", 
                     len(before), len(after), len(notfound), len(ignore))
        iofuncs.write_lines(
            before_file, sorted(zip(before, beforedates), key=lambda x:x[1]))
        iofuncs.write_lines(
            after_file, sorted(zip(after, afterdates), key=lambda x:x[1]))
        iofuncs.write_lines(notfound_file, notfound)
    return before, after, notfound


def compare_wang2007_query(test=False):
    """We compare MScanner retrieval to the complex PubMed query
    used in the Wang2007 paper."""
    dataset = "wang_2004"
    outdir = rc.working / "query" / dataset
    if not outdir.exists(): outdir.makedirs()
    if test:
        mindate = 19920101
        maxdate = 20020101
        t_mindate = 19920101
        t_maxdate = 20020101
        pos_file = rc.corpora / "Test" / "gdsmall.txt"
        neg_file = rc.articlelist
    else:
        mindate = 20040101
        maxdate = 20041231
        t_mindate = 20030101
        t_maxdate = 20031231
        pos_file = rc.corpora / "Wang2007" / "combined_pos.txt"
        neg_file = rc.corpora / "Wang2007" / "combined_neg.txt"
    all_pos = list(iofuncs.read_pmids(pos_file))
    all_neg = list(iofuncs.read_pmids(neg_file))
    env = Databases()
    old_pos, new_pos, nf_pos = split_by_date(env.artdb, all_pos, mindate, maxdate,
        outdir/"pos_before.txt", outdir/"pos_after.txt", outdir/"pos_notfound.txt")
    old_neg, new_neg, nf_neg = split_by_date(env.artdb, all_neg, mindate, maxdate,
        outdir/"neg_before.txt", outdir/"neg_after.txt", outdir/"neg_notfound.txt")
    limit = len(new_pos) + len(new_neg)
    QM = QueryManager(outdir, dataset, limit, None, env,
                      mindate, maxdate, t_mindate, t_maxdate)
    QM.query(old_pos, train_exclude=all_pos)
    new_results = [p for s,p in QM.results]
    truepos = result_relevance(new_results, set(new_pos))
    plot_rank_performance(outdir/"perf_vs_rank.png", 
                          outdir/"perf_vs_rank.txt",
                          truepos, len(new_pos),
                          compare_irrelevant=len(new_neg))
    QM.write_report(maxreport=1000)
    env.close()
    


if __name__ == "__main__":
    iofuncs.start_logger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression to execute"
    else:
        eval(sys.argv[1])
    logging.shutdown()
