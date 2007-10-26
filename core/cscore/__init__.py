"""Calculates citation scores"""

from __future__ import division
from mscanner.configuration import rc
import logging as log
import numpy as nx


                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


def pyscore(docs, featscores, offset, limit, threshold=None, exclude=[]):
    """Get scores for given documents
    
    We iterates over docs to yield scores. Skips members of exclude, and
    returns up to to limit results.
    
    @param docs: Iterator over (integer doc ID, array of feature ID) pairs

    @param featscores: Array of feature scores (mapping feature ID to score)
    
    @param offset: Arbitrary amount to add to citation score
    
    @param limit: Max number of results to return
    
    @param exclude: PMIDs to exclude from scoring
    
    @return: Iteration of (score, PMID) pairs
    """
    results = [(-100000, 0)] * limit
    import heapq
    ndocs = 0
    log.debug("Calculating article scores")
    marker = 0
    for idx, (docid, features) in enumerate(docs):
        if idx == marker:
            log.debug("Scored %d citations so far", idx)
            marker += 100000
        score = offset + nx.sum(featscores[features])
        if (threshold is None or score >= threshold) and docid not in exclude:
            ndocs += 1
            if score >= results[0][0]:
                heapq.heapreplace(results, (score,docid))
    if ndocs < limit:
        limit = ndocs
    return heapq.nlargest(limit, results)


def pyscore_adaptor(docstream, numdocs, featscores, offset,
                    limit, safety, threshold=None, exclude=[]):
    """Calls pyscore, given arguments suitable for cscore_pipe/dll"""
    from mscanner.medline.FeatureDatabase import FeatureStream
    docs = FeatureStream(open(docstream, "rb"))
    try:
        return pyscore(docs, featscores, offset, limit, threshold, exclude)
    finally:
        docs.close()


def cscore_pipe(docstream, numdocs, featscores, offset,
                limit, safety, threshold=None, exclude=[]):
    """Calculate article scores by piping to the cscore program
    
    The cscore program processes a feature stream to return a (score, pmid)
    pairs as a binary stream.
    
    @param docstream: Path to file containing feature vectors for documents to
    score, in L{mscanner.medline.FeatureDatabase.FeatureStream} format.
    
    @param numdocs: Number of documents in the stream of feature vectors.
    
    @param featscores: Vector of feature score doubles
    
    @param offset: Arbitrary amount to add to citation score

    @param limit: Maximum number of results to return.
    
    @param safety: Number of spare results in processing (because
    some might be members of exclude)
    
    @param exclude: PMIDs to remove from cscore results

    @param threshold: Cutoff score for including an article in the results

    @return: List of (score, PMID) pairs in decreasing order of score
    """
    import struct
    import subprocess as sp
    p = sp.Popen(
        [rc.cscore_path, docstream, str(numdocs), 
         str(len(featscores)), str(limit+safety), str(offset)],
        stdout=sp.PIPE, stdin=sp.PIPE)
    p.stdin.write(featscores.tostring())
    s = p.stdout.read(8)
    count = 0
    # Go through results in decreasing order to filter them
    while s != "":
        score, pmid = struct.unpack("fI", s)
        if (threshold is None or score > threshold) and pmid not in exclude:
            yield score, pmid
            count += 1
            if count >= limit:
                break
        s = p.stdout.read(8)
    p.stdout.close()


def cscore_dll(docstream, numdocs, featscores, offset,
               limit, safety, threshold=None, exclude=[]):
    """Calculate article scores, using ctypes to call cscores
    
    @param docstream: Path to file containing feature vectors for documents to
    score (formatted as in mscanner.medline.FeatureDatabase.FeatureStream)
    
    @param numdocs: Number of documents in the stream of feature vectors.
    
    @param featscores: Vector of feature score doubles
    
    @param offset: Arbitrary amount to add to citation score
    
    @param limit: Maximum number of results to return.
    
    @param safety: Number of spare results in processing (because
    some might be members of exclude)
    
    @param exclude: PMIDs to remove from cscore results

    @param threshold: Cutoff score for including an article in the results

    @return: List of (score, PMID) pairs in decreasing order of score
    """
    from ctypes import cdll, c_int, c_char_p, c_float, c_double
    from itertools import izip
    import numpy as nx
    # Set up arguments and call cscore2 function using ctypes
    carray = lambda x: nx.ctypeslib.ndpointer(dtype=x, ndim=1, flags='CONTIGUOUS')
    cscore = cdll.LoadLibrary(rc.cscore_dll)
    cscore.cscore.argtypes = [ 
        c_char_p, c_int, c_int, c_int, c_float,
        carray(nx.float64), carray(nx.float32), carray(nx.int32) ]
    o_scores = nx.zeros(limit+safety, dtype=nx.float32)
    o_pmids = nx.zeros(limit+safety, dtype=nx.int32)
    cscore.cscore(
        docstream, numdocs, len(featscores), limit+safety, offset,
        featscores, o_scores, o_pmids)
    # Go through results in decreasing order to filter them
    count = 0
    for score, pmid in izip(o_scores, o_pmids):
        if (threshold is None or score > threshold) and pmid not in exclude:
            yield score, pmid
            count += 1
            if count >= limit:
                break


score = pyscore_adaptor
"""Default score calculation function (parameters as for L{cscore_pipe})"""

def choose_score():
    """Select the fastest available score calculator and assign it
    to the module variable L{score}"""
    global score
    try:
        import ctypes
        if rc.cscore_dll.isfile():
            score = cscore_dll
    except ImportError:
        pass
    if score == pyscore_adaptor and rc.cscore_path.isfile():
        score = cscore_pipe

