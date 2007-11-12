"""Calculates the number of occurrences of each feature for
all articles in Medline between two dates."""

from __future__ import division
import logging as log
import numpy as nx
from path import path

from mscanner import update
from mscanner.configuration import rc
from mscanner.medline.FeatureDatabase import FeatureStream


class FeatureCounter:
    """Class for calculating feature counts in a subset of Medline"""
    
    counter_path = path(__file__).dirname() / "_FeatureCounter"
    
    def __init__(self,
                 docstream,
                 numdocs,
                 numfeats,
                 mindate=11110101,
                 maxdate=33330303,
                 exclude=set(),
                 ):
        update(self, locals())


    def py_counts(s):
        """Simply iterate over the documents and count how
        many times each feature occurs in the specified range"""
        featcounts = nx.zeros(s.numfeats, nx.int32)
        docs = FeatureStream(open(s.docstream, "rb"))
        try:
            for docid, date, features in docs:
                if (date >= s.mindate and date <= s.maxdate 
                    and docid not in s.exclude):
                    featcounts[features] += 1
        finally:
            docs.close()
        return featcounts


    def c_counts(s):
        """Pipes parameters to a C program that parses the stream
        of documents with features, which counts the number
        of occurrences of each feature, only considering documents
        added to Medline in the specified date range."""
        import subprocess as sp
        p = sp.Popen([
            s.counter_path, 
            s.docstream,
            str(s.numdocs),
            str(s.numfeats),
            str(s.mindate),
            str(s.maxdate),
            str(len(s.exclude)),
            ], stdout=sp.PIPE, stdin=sp.PIPE)
        p.stdin.write(nx.array(sorted(s.exclude)))
        return nx.fromfile(p.stdout, nx.int32, s.numfeats)
