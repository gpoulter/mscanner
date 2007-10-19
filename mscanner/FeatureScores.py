"""Calculates feature scores from occurrence counts"""

from __future__ import division
import logging as log
import numpy as nx

from mscanner import utils
from mscanner.Storage import Storage


                                     
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


class FeatureScores(object):
    """Feature score calculation and saving, with choice of calculation method,
    and methods to exclude certain kinds of features.
    
    @group featmap, pseudocount, mask, make_scores, get_postmask:
    Constructor parameters
    
    @ivar featmap: L{FeatureMapping} object

    @ivar pseudocount: Either a float which is the pseudocount for all
    features, or numpy array of per-feature counts, or None which triggers
    once-off creation of per-features counts equal to Medline frequency if
    Bayesian calculation is in use.
    
    @ivar mask: Boolean array for masking some features scores to zero
    (this is to exclude features by category, not by score).
    
    @ivar make_scores: Method used to calculate the feature scores.
    
    @ivar get_postmask: Method used to calculate a dynamic mask
    array once the feature scores are known.
    
    
    @group pos_counts, neg_counts, pdocs, ndocs: Set by L{update}
    
    @ivar pos_counts: Array of feature counts in positive documents
    
    @ivar neg_counts: Array of feature counts in negatives documents
    
    @ivar pdocs: Number of positive documents
    
    @ivar ndocs: Number of negative documents
    
    
    @group scores, pfreqs, nfreqs, offset: Set by L{make_scores}
    
    @ivar scores: Score of each feature
  
    @ivar pfreqs: Numerator of score fraction
    
    @ivar nfreqs: Denominator of score fraction
    
    @ivar offset: Value to be added to all article scores
   
    """

    def __init__(self, 
                 featmap,
                 pseudocount=None,
                 mask=None,
                 make_scores="scores_newpseudo",
                 get_postmask=None):
        """Initialise FeatureScores object (parameters are instance variables)"""
        self.featmap = featmap
        self.pseudocount = pseudocount
        self.mask = mask
        self.make_scores = getattr(self, make_scores)
        if isinstance(get_postmask, basestring):
            self.get_postmask = getattr(self, get_postmask)
        else:
            self.get_postmask = None


    def __len__(self):
        """Length is the number of features"""
        return len(self.featmap)


    def update(self, pos_counts, neg_counts, pdocs, ndocs):
        """Change the feature counts and numbers of documents, clear
        old score calculations, and calculate new scores."""
        utils.update(self, locals())
        self.offset = 0
        self.make_scores()
        self._mask_scores()
        for aname in ["_stats","_tfidf"]:
            try:
                delattr(self, aname)
            except AttributeError:
                pass


    def _mask_scores(self):
        """Set some feature scores to zero, effectively excluding them
        from consideration.  Uses L{mask} and L{get_postmask}"""
        if self.mask is not None:
            self.pfreqs[self.mask] = 0
            self.nfreqs[self.mask] = 0
            self.scores[self.mask] = 0
        if self.get_postmask:
            self.scores[self.get_postmask()] = 0


    def scores_noadjust(s):
        """Like L{scores_withabsence}, but neglects to adjust the feature scores.
        
        The motivation is that mysteriously we get a drop in averaged precision,
        so I'm seeing if the offset alone does anything, using the old feature scores."""
        s.scores_withabsence()
        s.pfreqs = s._pfreqs
        s.nfreqs = s._nfreqs
        s.scores = nx.log(s.pfreqs) - nx.log(s.nfreqs)
        print nx.sum(s.pfreqs), nx.sum(s.nfreqs)


    def scores_withabsence(s):
        """We model the document as a vector of occurring and non-occurring terms.
        
        The method pretends the probability of non-occurrence is 1, which is
        statistically incorrect in the multivariate Bernoulli model (the
        probability of a document with no features is 1 - adding features
        can only make the document more likely).
        """
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        s._pfreqs = (s.pseudocount+s.pos_counts) / (1+s.pdocs)
        s._nfreqs = (s.pseudocount+s.neg_counts) / (1+s.ndocs)
        s.pfreqs = s._pfreqs / (1-s._pfreqs)
        s.nfreqs = s._nfreqs / (1-s._nfreqs)
        s.scores = nx.log(s.pfreqs) - nx.log(s.nfreqs)
        s._constant_offset = nx.sum(nx.log(1-s._pfreqs)-nx.log(1-s._nfreqs))
        s._bayes_prior = nx.log(s.pdocs) - nx.log(s.ndocs)
        s.offset = s._constant_offset + s._bayes_prior


    def scores_newpseudo(s):
        """Use a pseudocount vector, or a constant pseudocount to get
        posterior feature probablilities.
        
        @return: Probability vectors for positive and negative features
        """
        s.offset = 0
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        s.pfreqs = (s.pseudocount+s.pos_counts) / (1+s.pdocs)
        s.nfreqs = (s.pseudocount+s.neg_counts) / (1+s.ndocs)
        s.scores = nx.log(s.pfreqs) - nx.log(s.nfreqs)


    def scores_oldpseudo(s):
        """Use a pseudocount vector, or a constant pseudocount to get
        posterior feature probablilities.   Instead of +1 in the
        denominator, we use +2*pseudocount.
        
        @deprecated: Performs poorly when number of inputs is small.
        
        @return: Probability arrays for positive and negative features"""
        s.offset = 0
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        s.pfreqs = (s.pseudocount+s.pos_counts) / (2*s.pseudocount+s.pdocs)
        s.nfreqs = (s.pseudocount+s.neg_counts) / (2*s.pseudocount+s.ndocs)
        s.scores = nx.log(s.pfreqs) - nx.log(s.nfreqs)


    def scores_rubin(s):
        """Uses Daniel Rubin's frequency smoothing heuristic, which
        replaces zero-frequency with 1e-8"""
        s.offset = 0
        s.pfreqs = s.pos_counts / float(s.pdocs)
        s.nfreqs = s.neg_counts / float(s.ndocs)
        s.pfreqs[s.pfreqs == 0.0] = 1e-8
        s.nfreqs[s.nfreqs == 0.0] = 1e-8
        s.scores = nx.log(s.pfreqs) - nx.log(s.nfreqs)


    def make_rare_positives(s):
        """Mask for positive-scoring features that do not occur in the positive set
        @return: Boolean array for masked out features
        """
        return (s.scores > 0) & (s.pos_counts == 0)


    def mask_nonpositives(s):
        """Mask for features not represented in the positives
        @return: Boolean array for masked out features
        """
        return s.pos_counts == 0


    @property 
    def stats(self):
        """A Storage instance with statistics about the features
        
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
        except AttributeError: 
            pass
        s = Storage()
        s.offset = self.offset
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
        s = self
        s.tfidf
        if isinstance(s.pseudocount, float):
            pseudocount = nx.zeros_like(s.scores) + s.pseudocount
        else:
            pseudocount = s.pseudocount
        for t, score in sorted(
            enumerate(s.scores), key=lambda x:x[1], reverse=True):
            if s.mask is not None and s.mask[t]:
                continue
            stream.write(
                u'%.3f,%d,%d,%.2e,%.2e,%.2e,%d,%.2f,%s,"%s"\n' % 
                (s.scores[t], s.pos_counts[t], s.neg_counts[t], 
                 s.pfreqs[t], s.nfreqs[t], pseudocount[t], t,
                 s.tfidf[t], s.featmap[t][1], s.featmap[t][0]))



def FeatureCounts(nfeats, featdb, docids):
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