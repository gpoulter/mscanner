"""Performs calculation of term and citation scores"""

from __future__ import division

                                     
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

import logging as log
import numpy as nx
from mscanner.support.utils import selfupdate


class FeatureInfo(object):
    """Class to calculate and store all information about the
    features, their distribution, and their scores
    
    @note: FeatureInfo models a table with columns (pos_counts, neg_counts,
    pfreqs, nfreqs, mask, scores), and links to a FeatureMapping table (id,
    name, type, count) for Medline occurrence counts.
    
    @note: FeatureInfo can be re-used (e.g. in cross validation) when
    the data changes - but must be recreated when parameters change.

    @ivar featmap: featuremap.FeatureMapping object
    
    @ivar pos_counts: Array of feature counts in positive documents. 
    
    @ivar neg_counts: Array of feature counts in negatives documents
    
    @ivar pdocs: Number of positive documents
    
    @ivar ndocs: Number of negative documents
    
    @ivar pseudocount: Either a float which is the pseudocount for all
    features, or numpy array of per-feature counts, or None which triggers
    once-off creation of per-features counts equal to Medline frequency if
    Bayesian calculation is in use.
    
    @ivar mask: Boolean array for masking some features scores to zero
    (this is to exclude features by category, not by score).
    
    @ivar frequency_method: Method for calculating feature 
    probabilities (constructor takes method name).
    
    @ivar post_masker: Method for calculating mask after feature
    scores are known (constructor takes method name or None)
  
    @ivar pfreqs: Probability of each feature given positive article
    
    @ivar nfreqs: Probability of each feature given negative article
    """

    def __init__(self, 
                 featmap, 
                 pos_counts=None, 
                 neg_counts=None, 
                 pdocs=None, 
                 ndocs=None,
                 pseudocount=None, 
                 mask=None,
                 frequency_method="getProbabilitiesBayes",
                 post_masker=None):
        """Initialise FeatureInfo object (parameters are instance variables)
        """
        selfupdate()
        self.getFrequencies = getattr(self, frequency_method)
        self.getPostMask = getattr(self, post_masker) if post_masker else None
        
    def __len__(self):
        """Length is the number of features"""
        return len(self.featmap)
    
    def _reset(self):
        """Reset the properties on the object because something has changed"""
        for aname in ["_scores","pfreqs","nfreqs","_stats"]:
            try:
                delattr(self, aname)
            except AttributeError: 
                pass
        
    def updateFeatureScores(self, pos_counts, neg_counts, pdocs, ndocs):
        """Change the feature counts and numbers of documents, clear
        old score calculations, and calculate new scores."""
        self._reset()
        selfupdate()
    
    @property
    def scores(self):
        """Log likelihood support score of each feature. 
        
        Side effect of calling stores the first time is that feature
        probabilities are also placed in the pfreqs and nfreqs attributes.
        """
        try:
            return self._scores
        except AttributeError: pass
        pfreqs, nfreqs = self.getFrequencies()
        scores = nx.log(pfreqs / nfreqs)
        if self.mask is not None:
            pfreqs[self.mask] = 0
            nfreqs[self.mask] = 0
            scores[self.mask] = 0
        self.pfreqs = pfreqs
        self.nfreqs = nfreqs
        self._scores = scores
        if self.getPostMask:
            self._scores[self.getPostMask()] = 0
        return self._scores

    def getProbabilitiesRubin(s):
        """Uses Daniel Rubin's frequency smoothing heuristic, which
        replaces zero-frequency with 1e-8"""
        pfreqs = s.pos_counts / float(s.pdocs)
        nfreqs = s.neg_counts / float(s.ndocs)
        pfreqs[pfreqs == 0.0] = 1e-8
        nfreqs[nfreqs == 0.0] = 1e-8
        return pfreqs, nfreqs

    def getProbabilitiesBayes(s):
        """Use a pseudocount vector, or a constant pseudocount to get
        posterior feature probablilities.
        """
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        pfreqs = (s.pos_counts+s.pseudocount) / (s.pdocs+1)
        nfreqs = (s.neg_counts+s.pseudocount) / (s.ndocs+1)
        return pfreqs, nfreqs

    def getProbabilitiesOldBayes(s):
        """Use a pseudocount vector, or a constant pseudocount to get
        posterior feature probablilities.   Instead of +1 in the
        denominator, we use +2*pseudocount.
        
        @deprecated: Performs poorly when number of inputs is small.
        """
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        pfreqs = (s.pos_counts+s.pseudocount) / (s.pdocs+2*s.pseudocount)
        nfreqs = (s.neg_counts+s.pseudocount) / (s.ndocs+2*s.pseudocount)
        return pfreqs, nfreqs

    def maskRarePositives(s):
        """Positive-scoring features not occurring in the positive set"""
        return (s.scores > 0) & (s.pos_counts == 0)
        
    def maskNonPositives(s):
        """Get rid of any features not represented in the positives"""
        return s.pos_counts == 0

    @property 
    def stats(self):
        """Storage object with statistics about the features
        
        Property has the following keys:
            - pos_occurrences: Total feature occurrences in positives
            - neg_occurrences: Total feature occurrences in negatives
            - feats_per_pos: Number of features per positive article
            - feats_per_neg: Number of features per negative article
            - distinct_feats: Number of distinct features
            - pos_distinct_feats: Number of of distinct features in positives
            - neg_distinct_feats: Number of of distinct features in negatives
        """
        try: 
            return self._stats
        except AttributeError: pass
        from mscanner.support.storage import Storage
        s = Storage()
        s.pdocs = self.pdocs
        s.ndocs = self.ndocs
        s.num_feats = len(self)
        s.pos_occurrences = int(nx.sum(self.pos_counts)) 
        s.feats_per_pos = 0.0
        if self.pdocs > 0:
            s.feats_per_pos = s.pos_occurrences / s.pdocs 
        s.neg_occurrences = int(nx.sum(self.neg_counts))
        s.feats_per_neg = 0.0
        if self.ndocs > 0:
            s.feats_per_neg = s.neg_occurrences / s.ndocs 
        s.pos_distinct_feats = len(nx.nonzero(self.pos_counts)[0]) 
        s.neg_distinct_feats = len(nx.nonzero(self.neg_counts)[0]) 
        self._stats = s
        return self._stats
        
    @property
    def tfidf(self):
        """TF-IDF scores for each feature
        
        Property calculates and caches TF-IDF scores for terms, where for term
        frequency (TF) we treat the positive corpus as a single large document,
        but for inverse document frequency (IDF) each citation is a separate
        document."""
        try: 
            return self._tfidf
        except AttributeError: 
            pass
        self._tfidf = nx.zeros(len(self.pos_counts), dtype=float)
        # Document frequency
        docfreq_t = self.pos_counts+self.neg_counts
        # Number of documents
        N = self.pdocs+self.ndocs # number of documents
        # Inverse Document Frequency (log N/df_t)
        u = (docfreq_t != 0)
        idf = nx.log(N / docfreq_t[u])
        # Term frequency
        tf = (self.pos_counts[u] / nx.sum(self.pos_counts))
        # Calculate TF.IDF
        self._tfidf[u] = tf * idf
        return self._tfidf

    def writeScoresCSV(self, stream):
        """Write features scores as CSV to an output stream"""
        stream.write(u"score,positives,negatives,numerator,"\
                     u"denominator,pseudocount,termid,tfidf,type,term\n")
        self.scores
        self.tfidf
        if isinstance(self.pseudocount, float):
            pseudocount = nx.zeros_like(self.scores) + self.pseudocount
        else:
            pseudocount = self.pseudocount
        s = self
        for t, score in sorted(
            enumerate(s.scores), key=lambda x:x[1], reverse=True):
            if s.mask is not None and s.mask[t]:
                continue
            stream.write(
                u'%.3f,%d,%d,%.2e,%.2e,%.2e,%d,%.2f,%s,"%s"\n' % 
                (s.scores[t], s.pos_counts[t], s.neg_counts[t], 
                 s.pfreqs[t], s.nfreqs[t], pseudocount[t], t,
                 s.tfidf[t], s.featmap[t][1], s.featmap[t][0]))



def countFeatures(nfeats, featdb, docids):
    """Count occurrenes of each feature in a set of articles

    @param nfeats: Number of distinct features (size of array)

    @param featdb: Mapping from document ID to list of feature IDs

    @param docids: List of document IDs whose features are to be counted

    @return: Array of length nfeats, with occurrence count of each feature
    """
    counts = nx.zeros(nfeats, nx.int32)
    for docid in docids:
        counts[featdb[docid]] += 1
    return counts


def iterScores(docs, featscores, limit, threshold=None, exclude=[]):
    """Iterate over docs, yielding scores (calculated using featscores array),
    while skipping over members of exclude, returning up to limit results.
    
    @param docs: Iterator over (integer doc ID, array of feature ID) pairs

    @param featscores: Array of feature scores (mapping feature ID to score)
    
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
        score = nx.sum(featscores[features])
        if (threshold is None or score >= threshold) and docid not in exclude:
            ndocs += 1
            if score >= results[0][0]:
                heapq.heapreplace(results, (score,docid))
    if ndocs < limit:
        limit = ndocs
    return heapq.nlargest(limit, results)


def iterCScores(progpath, docstream, numdocs, featscores, 
                limit, safety, threshold=None, exclude=[]):
    """The cscore program processes the feature stream to return a binary
    string of (score, pmid) pairs. This function calls it and converts the
    output to python objects, and filters occurrences of input citations.
    
    @param progpath: Path to the cscore program
    
    @param docstream: Path to file containing feature vectors for documents to
    score (formatted as in mscanner.featuredb.FeatureStream)
    
    @param numdocs: Number of documents in the stream of feature vectors.
    
    @param featscores: Vector of feature score doubles
    
    @param limit: Maximum number of results to return.
    
    @param safety: Number of additional results the cscore should return (since
    some might be members of exclude)
    
    @param exclude: PMIDs to remove from cscore results

    @param threshold: Cutoff score for including an article in the results

    @return: List of (score, PMID) pairs in decreasing order of score
    """
    import struct
    import subprocess as sp
    p = sp.Popen(
        [progpath, docstream, str(numdocs), 
         str(len(featscores)), str(limit+safety)],
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


def iterCScores2(docstream, numdocs, featscores, 
                 limit, safety, threshold=None, exclude=[]):
    """Use ctypes to access cscores for calculating feature scores.
    
    @param docstream: Path to file containing feature vectors for documents to
    score (formatted as in mscanner.featuredb.FeatureStream)
    
    @param numdocs: Number of documents in the stream of feature vectors.
    
    @param featscores: Vector of feature score doubles
    
    @param limit: Maximum number of results to return.
    
    @param safety: Number of additional results the cscore should return (since
    some might be members of exclude)
    
    @param exclude: PMIDs to remove from cscore results

    @param threshold: Cutoff score for including an article in the results

    @return: List of (score, PMID) pairs in decreasing order of score
    """
    from ctypes import cdll, c_int, c_char_p
    from itertools import izip
    import numpy as nx
    # Set up arguments and call cscore2 function using ctypes
    carray = lambda x: nx.ctypeslib.ndpointer(dtype=x, ndim=1, flags='CONTIGUOUS')
    cscore = cdll.LoadLibrary(r"cscore2.dll")
    cscore.cscore.argtypes = [ 
        c_char_p, c_int, c_int, c_int, carray(nx.float64),
        carray(nx.float32), carray(nx.int32) ]
    o_scores = nx.zeros(limit+safety, dtype=nx.float32)
    o_pmids = nx.zeros(limit+safety, dtype=nx.int32)
    cscore.cscore(
        docstream, numdocs, len(featscores), limit+safety, 
        featscores, o_scores, o_pmids)
    # Go through results in decreasing order to filter them
    count = 0
    for score, pmid in izip(o_scores, o_pmids):
        if (threshold is None or score > threshold) and pmid not in exclude:
            yield score, pmid
            count += 1
            if count >= limit:
                break


def retrievalTest(results, test):
    """Perform a retrieval test of query results against a gold standard.

    The false positives can be calculated as rank minus true positives.

    @param results: List of result PubMed IDs
    
    @param test: Set of gold-standard articles to look for in results

    @return: True positives as function of rank (measured with respect
    to the test set)
    """
    assert isinstance(test, set)
    TP_total = nx.zeros(len(results), nx.int32)
    for idx, pmid in enumerate(results):
        TP_total[idx] = TP_total[idx-1] if idx > 0 else 0
        if pmid in test:
            TP_total[idx] += 1
    return TP_total
