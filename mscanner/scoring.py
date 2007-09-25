"""Performs calculation of term scores"""

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
from mscanner.support.utils import update


class FeatureInfo(object):
    """Class to calculate and store all information about the
    features, their distribution, and their scores
    
    @note: FeatureInfo is like a table with columns of pos_counts, neg_counts,
    pfreqs, nfreqs, mask, scores.  It links to a L{FeatureMapping} table (id,
    name, type, count) for Medline occurrence counts.
    
    FeatureInfo can be re-used (e.g. in cross validation) when the data changes
    but must be recreated when rc parameters change.

    @ivar featmap: L{FeatureMapping} object
    
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
        """Initialise FeatureInfo object (parameters are instance variables)"""
        update(self, locals())
        self.getFrequencies = getattr(self, frequency_method)
        self.getPostMask = getattr(self, post_masker) if post_masker else None


    def __len__(self):
        """Length is the number of features"""
        return len(self.featmap)


    def _reset(self):
        """Reset properties because something changed"""
        for aname in ["_scores","pfreqs","nfreqs","_stats","_tfidf"]:
            try:
                delattr(self, aname)
            except AttributeError: 
                pass


    def update_features(self, pos_counts, neg_counts, pdocs, ndocs):
        """Change the feature counts and numbers of documents, clear
        old score calculations, and calculate new scores."""
        self._reset()
        update(self, locals())


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
        
        @return: Probability vectors for positive and negative features
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

        @return: Probability arrays for positive and negative features"""
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        pfreqs = (s.pos_counts+s.pseudocount) / (s.pdocs+2*s.pseudocount)
        nfreqs = (s.neg_counts+s.pseudocount) / (s.ndocs+2*s.pseudocount)
        return pfreqs, nfreqs


    def maskRarePositives(s):
        """Mask for positive-scoring features that do not occur in the positive set
        @return: Boolean array for masked out features
        """
        return (s.scores > 0) & (s.pos_counts == 0)


    def maskNonPositives(s):
        """Mask for features not represented in the positives
        @return: Boolean array for masked out features
        """
        return s.pos_counts == 0


    @property 
    def stats(self):
        """A storage instance with statistics about the features
        
        The following keys are present:
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
        """Vector of TF-IDF scores for each feature
        
        Cache TF-IDF scores for terms, where for term frequency (TF) we treat
        the positive corpus as a single large document, but for inverse
        document frequency (IDF) each citation is a separate document."""
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


    def get_best_tfidfs(self, count):
        """Construct a table about the terms with the bets TF.IDF
        
        @param count: Number of rows to return
        
        @return: List of 
        (Term ID, TFIDF, (term, term_type), term score, pos count, neg count)
        """
        from heapq import nlargest
        best_tfidfs = nlargest(
            count, enumerate(self.tfidf), key=lambda x:x[1])
        return [ (t, tfidf, self.featmap[t], self.scores[t], 
                  self.pos_counts[t], self.neg_counts[t])
                  for t, tfidf in best_tfidfs ]


    def write_csv(self, stream):
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



def count_features(nfeats, featdb, docids):
    """Count occurrenes of each feature in a set of articles

    @param nfeats: Number of distinct features (length of L{docids})

    @param featdb: Mapping from document ID to array of feature IDs

    @param docids: Iterable of document IDs whose features are to be counted

    @return: Array of length L{nfeats}, containing occurrence count of each feature
    """
    counts = nx.zeros(nfeats, nx.int32)
    for docid in docids:
        counts[featdb[docid]] += 1
    return counts