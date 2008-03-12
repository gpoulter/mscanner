"""Calculates feature scores from occurrence counts"""

from __future__ import division
import logging
import numpy as nx
from itertools import izip
from heapq import nlargest

from mscanner import update, delattrs
from mscanner.core.Storage import Storage


                                     
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
    
    @group Set via constructor: featmap, scoremethod
    
    @ivar featmap: L{FeatureMapping} object
    
    @ivar scoremethod: Name of the method to call for calculating feature scores.
    
    
    @group Set by update: pdocs, ndocs, pos_counts, neg_counts, selected,
    features, pos_selected, neg_selected, prior
    
    @ivar pdocs: Number of positive documents.
    
    @ivar ndocs: Number of negative documents.
    
    @ivar pos_counts: Array with number of positive occurrences for each feature.
    
    @ivar neg_counts: Array with number of negative occurrences for each feature.
    
    @ivar selected: Boolean array with True position corresponds to selected features.
    
    @ivar features: Array listing the feature IDs of the L{selected} features.

    @ivar pos_selected: Array with number of positive documents for each
    selected feature (corresponds with L{features}).
    
    @ivar neg_selected: Array with number of negative documents for each
    selected feature (corresponds with L{features}).

    @ivar prior: The prior log odds of relevance.  If None, estimate
    using the log of the ratio of relevant to irrelevant articles in the data.


    @group Set via scoremethod: scores, pfreqs, nfreqs, base
    
    @ivar scores: Score of each feature.
  
    @ivar pfreqs: For positive documents, the estimated fraction having each feature.
    
    @ivar nfreqs: For negative documents, the estimated fraction having each feature.
    
    @ivar base: Value to be added to all article scores (typically the
    score of an article with zero features).
    """

    def __init__(self, featmap, scoremethod):
        """Initialise FeatureScores object (parameters are instance variables)"""
        update(self, locals())


    @staticmethod
    def Defaults(featmap):
        """Create a FeatureScores parameterised using the RC defaults.
        @param featmap: L{FeatureMapping} instance to use."""
        from mscanner.configuration import rc
        return FeatureScores(featmap, rc.scoremethod)


    def scores_of(self, docs):
        """Calculate scores of a list of documents.
        @param docs: List/iterable of feature vectors.
        @return: Vector of scores corresponding to L{docs}.
        """
        return self.base + self.prior +\
               nx.array([nx.sum(self.scores[d]) for d in docs], nx.float32)


    def update(self, pos_counts, neg_counts, pdocs, ndocs, prior=None):
        """Change the feature counts and numbers of documents and re-calculate
        the feature scores."""
        # Set attributes from parameters
        self.pdocs = pdocs
        self.ndocs = ndocs
        assert len(pos_counts) == len(neg_counts)
        self.pos_counts = pos_counts
        self.neg_counts = neg_counts
        self.prior = prior
        if self.prior is None:
            if pdocs == 0 or ndocs == 0:
                self.prior = 0.0
            else:
                # Convert prior probability of relevance to a log odds
                self.prior = nx.log(pdocs/ndocs)
        # Get boolean array for selected features
        self.selected = self.get_selected()
        # Get list of unmasked features and condensed count vectors
        self.features = nx.arange(len(self.selected))[self.selected]
        self.pos_selected = self.pos_counts[self.selected]
        self.neg_selected = self.neg_counts[self.selected]
        # Call the method to calculate feature scores (it respects the mask)
        getattr(self, self.scoremethod)()
        # Get rid of TFIDF property to recalculate
        delattrs(self, "_stats", "_tfidf")


    def get_selected(s):
        """Return a boolean feature selection array. True positions correspond
        to features that we keep.  Especially with word features, the
        number of selected features must be kept reasonably small.
        """
        from mscanner.configuration import rc
        selected = nx.ones(len(s.pos_counts), nx.bool)
        if rc.type_mask:
            # Keep only features that are not of the specified types
            selected &= ~s.featmap.type_mask(rc.type_mask)
        if rc.mincount > 0:
            # Keep only features having enough occurrences
            selected &= ((s.pos_counts + s.neg_counts) >= rc.mincount)
        if rc.positives_only:
            # Keep only features found in a relevant article
            selected &= (s.pos_counts != 0)
        if rc.min_infogain > 0:
            # Keep features by information gain. Infogain is expensive to 
            # calculate, so do it only for features that have passed so far.
            selected[selected] &= (s.infogain(selected) >= rc.min_infogain)
        return selected


    def infogain(s, selected):
        """Calculate the Information Gain (L{infogain}) of so-far-selecte features.
        
        @param selected: A boolean array over all features, where True marks
        features that have passed feature selection so far.
        
        @return: Array over selected features with information gain values.
        """
        def S(p): return -p*nx.log2(p) # Entropy in bits
        logging.debug("Calculating info gain on %d features out of %d", sum(selected), len(selected))
        R = s.pdocs # relevant
        I = s.ndocs # irrelevant
        N = R+I # total 
        R1 = s.pos_counts[selected] # relevant and term present
        I1 = s.neg_counts[selected] # irrelevant and term present
        N1 = R1+I1 # term present
        
        # MAP est prob of term being present
        p1 = (N1+1)/(N+2)  
        # p0 = (N0+1)/(N+2) = 1-p1 with N0+N1=N
        
        # Entropy non-partitioned distribution of classes
        SC = S(R/N) + S(I/N)
        
        # Entropy of MAP est prob of relevance/irrelevance given term present
        SCg1 = S((R1+1)/(N1+2)) + S((I1+1)/(N1+2))
        
        # Entropy of MAP est prob of relevance/irrelevance given term absent
        # SCg0 = S((R0+1)/(N0+2)) + S((I0+1)/(N0+2)) with R0+R1=R, I0+I1=I, N0+N1=N
        SCg0 = S((R+1-R1)/(N+2-N1)) + S((I+1-I1)/(N+2-N1)) 
        
        # Entropy diff between original and partitioned distribution of classes
        return SC - (p1*SCg1 + (1-p1)*SCg0)  # p0+p1=1
    


    def scores_bayes(s, pos_a, pos_ab, neg_a, neg_ab):
        """Estimate support scores of features assuming documents are generated
        by a multivariate Bernoulli distribution. For non-occurring features,
        we use a base score and adjust the score for feature occurrence.
        Evaluation of frequencies is done only on selected features (the
        feature ID for each item is the corresponding item in L{features}), but
        the final L{scores} array is full-length (feature ID of each item is
        the array index).
        
        @param pos_a, pos_ab: Beta prior (a=successes, ab=total) for relevant articles.
        @param neg_a, neg_ab: Beta prior (a=successes, ab=total) for irrelevant articles.
        """
        # Posterior probability of term given relevance of document
        s.pfreqs = (pos_a+s.pos_selected) / (pos_ab+s.pdocs) 
        # Posterior probability of term given irrelevance of document
        s.nfreqs = (neg_a+s.neg_selected) / (neg_ab+s.ndocs)
        # Support scores for bernoulli successes (then apply mask)
        s.success_scores = nx.log(s.pfreqs / s.nfreqs)
        # Support scores for bernoulli failures (then apply mask)
        s.failure_scores = nx.log((1-s.pfreqs) / (1-s.nfreqs))
        # Convert success/failure to base score and occurrence score
        s.base = nx.sum(s.failure_scores)
        # Create score vector for ALL features, single-precision, with zeros at non-selected features
        s.scores = nx.zeros(len(s.featmap), nx.float32)
        s.scores[s.selected] = nx.float32(s.success_scores - s.failure_scores)


    def scores_laplace(s):
        """For feature probabilities we use the Laplace prior of 1 success and
        1 failure in each class. Use with a mask that excludes features not
        occurring in relevant articles."""
        s.scores_bayes(1, 2, 1, 2)


    def scores_laplace_split(s):
        """For feature probabilities we use a Laplace prior, of 1 success and 1
        failure in total, split between the classes according to size. This
        avoids class skew problems."""
        p = s.pdocs / (s.pdocs + s.ndocs)
        s.scores_bayes(p, 2*p, 1-p, 2*(1-p))


    def scores_bgfreq(s):
        """The prior is 'background frequency' successes, out of 1 total occurrence 
        in each class."""
        bgvec = s.featmap.counts[s.selected] / s.numdocs
        s.scores_bayes(bgvec, 1, bgvec, 1)


    def scores_rubin(s):
        """Maximum likelihood classifier: score is a product of the likelihood
        ratios for each occurring features. Uses ML estimates, with zero
        probabilities are replaced by 1e-8 (ignores other masks)"""
        s.base = 0
        s.prior = 0
        s.pfreqs = s.pos_counts / float(s.pdocs)
        s.nfreqs = s.neg_counts / float(s.ndocs)
        s.pfreqs[s.pfreqs == 0.0] = 1e-8
        s.nfreqs[s.nfreqs == 0.0] = 1e-8
        s.scores = nx.log(s.pfreqs) - nx.log(s.nfreqs)
        s.scores[~s.selected] = 0


    @property 
    def stats(self):
        try: 
            return self._stats
        except AttributeError: 
            self._stats = self.FeatureStats(self)
            return self._stats


    class FeatureStats:
        """Stores statistics about a L{FeatureScores}
        
        @ivar feats_selected: Number of features that passed feature selection.
        @ivar feats_total: Total number of features in the feature space.
        
        @ivar pos_docs, neg_docs: Number of documents in positive/negative
        corpus.
        
        @ivar pos_occurrences, neg_occurrences: Total number of feature
        occurrences in the positive/negative corpus.
        
        @ivar pos_average, neg_average: Average number of features per
        positive/negative document.
        
        @ivar pos_distinct, neg_distinct: Number of different unmasked features
        in positive/negative document.
        """
        
        def __init__(s, featscores):
            fs = featscores
            s.feats_selected = sum(fs.selected)
            s.feats_total = len(fs.selected)
            s.pos_docs = fs.pdocs
            s.neg_docs = fs.ndocs
            s.pos_occurrences = int(nx.sum(fs.pos_selected))
            s.neg_occurrences = int(nx.sum(fs.neg_selected))
            s.pos_average = s.pos_occurrences / fs.pdocs if fs.pdocs > 0 else 0.0
            s.neg_average = s.neg_occurrences / fs.ndocs if fs.ndocs > 0 else 0.0
            s.pos_distinct = sum(fs.pos_selected != 0)
            s.neg_distinct = sum(fs.neg_selected != 0)


    @property
    def tfidf(self):
        """Vector of TF-IDF scores for each selected feature. To obtain term
        frequency (TF) we treat the positive corpus as a single large document.
        To obtain inverse document frequency (IDF), we consider each Medline
        record to be a separate document."""
        try: 
            return self._tfidf
        except AttributeError: 
            pass
        self._tfidf = nx.zeros(len(self.pos_selected), dtype=float)
        # Document frequency
        df_t = self.pos_selected + self.neg_selected
        # Number of documents
        N = self.pdocs + self.ndocs
        # Inverse Document Frequency (log N/df_t)
        idf = nx.log(N / df_t)
        # Get rid of zero-frequency features
        idf[df_t == 0] = 0
        # Term frequency
        tf = self.pos_selected / nx.sum(self.pos_selected)
        # Calculate TF.IDF
        self._tfidf = tf * idf
        return self._tfidf
    

    def get_best_tfidfs(s, maxfeats, min_tfidf=0.001):
        """Get a table of features with the highest TFIDF values, up to a limit
        and only for TF-IDFs above a minimum.
        
        @param maxfeats: Maximum number of features to return.
        
        @param min_tfidf: Minimum value of TF*IDF to return.
        
        @return: List of (Feature ID, TFIDF, (feature, featureclass), 
        feature score, positive count, negative count)."""
        best_tfidfs = nlargest(maxfeats, izip(s.tfidf, s.features))
        return [ (t, tfidf, s.featmap.get_feature(t), s.scores[t], 
                  s.pos_counts[t], s.neg_counts[t])
                  for tfidf, t in best_tfidfs if tfidf >= min_tfidf ]


    def write_csv(s, stream, maxfeats=None):
        """Write scores of selected features as CSV to an output stream.
        
        @param stream: File-like object, supporting a unicode write method.
        
        @param maxfeats: Write only the top N selected features (None for all selected).
        """
        logging.info("There are %d selected features out of %d total.", sum(s.selected), len(s.selected))
        # Output features by decreasing score
        stream.write(u"score,pos_count,neg_count,termid,type,term\n")
        if maxfeats is None: maxfeats = len(s.selected)
        for (score, t) in nlargest(maxfeats, izip(s.scores[s.selected], s.features)):
            fname, ftype = s.featmap.get_feature(t)
            stream.write(u'%.3f,%d,%d,%d,%s,"%s"\n' % 
            (s.scores[t], s.pos_counts[t], s.neg_counts[t], t, ftype, fname))
