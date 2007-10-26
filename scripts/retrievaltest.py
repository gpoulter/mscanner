#!/usr/bin/env python

"""Performs retrieval test analysis, where a subset of the input is used to
query, and the results are compared against the remainder of the input.

Usage::

    python retrievaltest.py <dataset> <filename>
    
For example::

    python retrievaltest.py pg07 'pharmgkb-070205.txt'
"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import random
import sys
from path import path

from mscanner.configuration import rc, start_logger
from mscanner.core.QueryManager import QueryManager
from mscanner.core.Plotter import Plotter


### RETRIEVAL TEST
## Name of file with list of testing PMIDs
rc.report_retrieval_test_pmids = path("retrieval_test.txt")
## Name of file with cumulative count of retrieved test PMIDs
rc.report_retrieval_stats = path("retrieval_stats.txt")
## Name of retrieval vs rank graph
rc.report_retrieval_graph = path("retrieval.png")


def compare_results_to_standard(results, test):
    """Test of query results against a gold standard.

    @param results: List of result PubMed IDs
    
    @param test: Set of gold-standard articles to look for in results

    @return: True positives as function of rank (measured with respect to the
    test set). False positives are equal to rank minus true positives."""
    assert isinstance(test, set)
    import numpy as nx
    TP_total = nx.zeros(len(results), nx.int32)
    for idx, pmid in enumerate(results):
        TP_total[idx] = TP_total[idx-1] if idx > 0 else 0
        if pmid in test:
            TP_total[idx] += 1
    return TP_total


class RetrievalTest(QueryManager):
    """Performs retrieval-testing analysis, in which a subset of the data is
    used to query, and the results are tested against the rest of the data."""
    
    def retrieval_query(self, input):
        """Splits the input into training and testing, plots how many testing
        PMIDs the query returns when trained on the training PMIDs.
        
        @return: Array with cumulative test PMIDs as function of rank
        """
        rc.threshold = None # No threshold for this stuff
        # Split the input into train and test sections
        self._load_input(input)
        if len(self.pmids) == 0: return
        pmids_list = list(self.pmids)
        random.shuffle(pmids_list)
        subset_size = int(rc.retrieval_test_prop * len(pmids_list))
        self.pmids = set(pmids_list[:subset_size])
        test_pmids = set(pmids_list[subset_size:])
        # Write test PMIDs
        (self.outdir/rc.report_retrieval_test_pmids).write_lines(
            [str(x) for x in test_pmids])
        # Get feature info and result scores
        self._make_feature_info()
        self._make_results()
        self._save_results()
        # Test the results against the test PMIDs
        cumulative = compare_results_to_standard(
            [p for s,p in self.results], test_pmids)
        # Write the number of TP at each rank
        (self.outdir/rc.report_retrieval_stats).write_lines(
            [str(x) for x in cumulative])
        # Graph TP vs FP (FP = rank-TP)
        plotter = Plotter()
        plotter.plot_retrieved_positives(
            self.outdir/rc.report_retrieval_graph, 
            cumulative, len(test_pmids))


def retrieval(dataset, filename):
    """Carry out retrieval test
    
    @param dataset: Name of the task
    @param filename: Name of corpus to load
    """
    rc.limit = 1000 # Force 1000 results
    rc.retrieval_test_prop = 0.2
    rc.dataset = dataset
    op = RetrievalTest(rc.web_report_dir / dataset)
    op.retrieval_query(rc.corpora / filename)
    op.env.close()


if __name__ == "__main__":
    start_logger()
    if len(sys.argv) != 3:
        print "Please provide dataset and input file name"
    else:
        retrieval(sys.argv[1], sys.argv[2])
