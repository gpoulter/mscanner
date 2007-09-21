#!/usr/bin/env python

"""Performs queries using the datasets from the MScanner paper

Choose operation by executing a Python expression::

    python query.py 'query("pg04","pg07")'
"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import random
import sys
from mscanner.configuration import rc, initLogger
from mscanner import plotting, scorefile, scoring, queryenv


class RetrievalTest(queryenv.Query):
    """Performs retrieval-testing analysis, in which a subset of the data is
    used to query, and the results are tested against the rest of the data."""
    
    def testRetrieval(self, input):
        """Splits the input into training and testing, plots how many testing
        PMIDs the query returns when trained on the training PMIDs.
        
        @return: Array with cumulative test PMIDs as function of rank
        """
        rc.threshold = None # No threshold for this stuff
        # Split the input into train and test sections
        self.load_pmids(input)
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
        self.make_featinfo()
        self.make_results()
        self.save_results()
        # Test the results against the test PMIDs
        cumulative = scoring.compare_results_to_standard(
            [p for s,p in self.results], test_pmids)
        # Write the number of TP at each rank
        (self.outdir/rc.report_retrieval_stats).write_lines(
            [str(x) for x in cumulative])
        # Graph TP vs FP (FP = rank-TP)
        plotter = plotting.Plotter()
        plotter.plotRetrievalGraph(
            self.outdir/rc.report_retrieval_graph, 
            cumulative, len(test_pmids))


def retrieval(dataset, filename):
    """Carry out retrieval test
    
    @param dataset: Name of the task
    @param filename: Name of corpus to load"""
    rc.limit = 1000 # Force 1000 results
    rc.retrieval_test_prop = 0.2
    rc.dataset = dataset
    op = RetrievalTest(rc.report_dir / dataset)
    op.testRetrieval(rc.corpora / filename)
    op.env.close()


if __name__ == "__main__":
    initLogger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression"
    else:
        eval(sys.argv[1])
