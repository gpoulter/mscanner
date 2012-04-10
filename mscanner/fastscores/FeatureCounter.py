"""Calculates the number of occurrences of each feature for
all articles in Medline between two dates."""

from __future__ import division
import numpy as nx
from path import path
import struct

from mscanner import update
from mscanner.medline.FeatureStream import FeatureStream


class FeatureCounter:
    """Class for calculating feature counts in a subset of Medline.

    @ivar docstream: Path to file containing feature vectors for documents to
    score, in FeatureStream format.
    
    @ivar numdocs: Number of documents in the stream article stream.
    
    @ivar numfeats: Number of distinct features in Medline (length of the 
    vector of feature counts).

    @ivar mindate: YYYYMMDD as integer: ignore documents before this date.
    
    @ivar maxdate: YYYYMMDD as integer: ignore documents after this date.

    @ivar exclude: PubMed IDs that are not allowed to appear in the results.
    """


    def __init__(self,
                 docstream,
                 numdocs,
                 numfeats,
                 mindate=None,
                 maxdate=None,
                 exclude=set(),
                 ):
        if mindate is None: mindate = 11110101
        if maxdate is None: maxdate = 99990101
        self.counter_path = path(__file__).dirname() / "_FeatureCounter"
        update(self, locals())


    def py_counts(s):
        """Simply iterate over the documents and count how
        many times each feature occurs in the specified range
        
        @return: Number of documents counted, and vector of feature counts."""
        featcounts = nx.zeros(s.numfeats, nx.int32)
        docs = FeatureStream(s.docstream, rdonly=True)
        ndocs = 0
        try:
            for docid, date, features in docs.iteritems():
                if (date >= s.mindate and date <= s.maxdate 
                    and docid not in s.exclude):
                    featcounts[features] += 1
                    ndocs += 1
        finally:
            docs.close()
        return ndocs, featcounts


    def c_counts(s):
        """Pipes parameters to a C program that parses the stream
        of documents with features, which counts the number
        of occurrences of each feature, only considering documents
        added to Medline in the specified date range.
        
        @return: Number of documents counted, and vector of feature counts."""
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
        # First integer of output is the number of documents parsed
        ndocs = struct.unpack("I", p.stdout.read(4))[0]
        # Then a vector of feature counts
        featcounts = nx.fromfile(p.stdout, nx.int32, s.numfeats)
        # Make child terminate
        p.stdout.close()
        return ndocs, featcounts
