"""Calculates feature scores from occurrence counts"""

from __future__ import division
import logging
import numpy as nx

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
    
    
    @group Set by update: pos_counts, neg_counts, pdocs, ndocs, prior
    
    @ivar pos_counts: Number of positive documents having each feature.
    
    @ivar neg_counts: Number of negative documents having each feature.
    
    @ivar pdocs: Number of positive documents.
    
    @ivar ndocs: Number of negative documents.
    
    @ivar prior: The prior log odds of relevance.  If None, estimate
    using the log of the ratio of relevant to irrelevant articles in the data.
    
    @ivar mask: Boolean array of features to exclude (for feature selection).
    Features corresponding to True positions are set to zero.


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


    def scores_of(self, featdb, pmids):
        """Calculate vector of scores given an iterable of PubMed IDs.
        
        @param featdb: Mapping from PMID to feature vector
        @param pmids: Iterable of keys into L{featdb}
        @return: Vector containing document scores corresponding to the pmids.
        """
        off = self.base + self.prior
        sc = self.scores
        return nx.array([off+nx.sum(sc[featdb[d]]) for d in pmids], nx.float32)


    def __len__(self):
        """Number of features"""
        return len(self.featmap)


    def update(self, pos_counts, neg_counts, pdocs, ndocs, prior=None):
        """Change the feature counts and numbers of documents and re-calculate
        the feature scores."""
        if prior is None:
            if pdocs == 0 or ndocs == 0:
                prior = 0.0
            else:
                # Convert prior probability of relevance to a log odds
                prior = nx.log(pdocs/ndocs)
        update(self, locals())
        # Call the method to exclude certain features
        self.mask = self.get_mask()
        # Call the method to calculate feature scores
        getattr(self, self.scoremethod)()
        # Allow TFIDF property to recalculate
        delattrs(self, "_stats", "_tfidf")


    def get_mask(self):
        """Create feature mask for excluded feature classes, document
        frequency cutoff and information gain cutoff - using RC parameters."""
        from mscanner.configuration import rc
        r = nx.zeros(len(self.pos_counts), nx.bool)
        if rc.class_mask:
            # Exclude certain classes of features (word/mesh/issn)
            r |= self.featmap.class_mask(rc.class_mask)
        if rc.mincount > 0:
            # Exclude features not having enough occurrences
            r |= ((self.pos_counts + self.neg_counts) < rc.mincount)
        if rc.positives_only:
            # Exclude features not occurring in any relevant articles
            r |= (self.pos_counts == 0)
        if rc.min_infogain > 0:
            # Exclude features not having enough information gain
            r |= (self.infogain() < rc.min_infogain)
        return r


    def infogain(s):
        """Calculate the Information Gain (L{infogain}) of each feature."""
        # Entropy in bits
        def info(p): 
            return -p*nx.log2(p)
        
        # Number of documents
        N = s.pdocs + s.ndocs 
        # Number documents that have term i 
        T = s.pos_counts + s.neg_counts 
        
        # Fraction of documents that have term i 
        # (probability of term being present)
        pT = T / N
        
        # Fraction of documents that are relevant
        # (prior probability of relevance)
        pR = s.pdocs / N
        
        # Fraction of relevant documents amonst documents that have term i
        # (posterior probability of relevance given presence of term)
        pRgT = (s.pos_counts+1) / (T+2)
        
        # Fraction of relevant documents amongst documents not having term i
        # (posterior probability of relevance given absence of term)
        pRgNT = (s.pdocs-s.pos_counts+1) / (N-T+2)
        
        pI = s.ndocs / N
        pIgT = (s.neg_counts+1) / (T+2)
        pIgNT = (s.ndocs-s.neg_counts+1) / (N-T+2)
        
        # Information gain is entropy before knowing whether it's present or
        # absent, minus expected entropy after knowing.
        return (info(pR)+info(pI)) - (
            pT*(info(pRgT)+info(pIgT)) + (1-pT)*(info(pRgNT)+info(pIgNT)))


    def scores_bayes(s, pos_a, pos_ab, neg_a, neg_ab):
        """Estimate support scores of features assuming documents are generated
        by a multivariate Bernoulli distribution. Applies the L{mask}
        attribute. Feature non-occurrence is modeled as a base score for the
        document with no features, and an adjustment to the feature occurrence
        scores.
        
        @param pos_a, pos_ab: Beta prior (a=successes, ab=total) for relevant articles.
        @param neg_a, neg_ab: Beta prior (a=successes, ab=total) for irrelevant articles.
        """
        # Posterior probability of term given relevance of document
        s.pfreqs = (pos_a+s.pos_counts) / (pos_ab+s.pdocs) 
        # Posterior probability of term given irrelevance of document
        s.nfreqs = (neg_a+s.neg_counts) / (neg_ab+s.ndocs)
        # Support scores for bernoulli successes (apply mask)
        s.success_scores = nx.log(s.pfreqs / s.nfreqs)
        s.success_scores[s.mask] = 0
        # Support scores for bernoulli failures (apply mask)
        s.failure_scores = nx.log((1-s.pfreqs) / (1-s.nfreqs))
        s.failure_scores[s.mask] = 0
        # Conver success/failure to base score and occurrence score
        s.base = nx.sum(s.failure_scores)
        s.scores = s.success_scores - s.failure_scores


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
        bgvec = s.background_vector
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


    @property
    def background_vector(s):
        """Vector of background probability of occurrence
        of each feature in Medline (for use in regularisation)."""
        if not hasattr(s, "_background_vector"):
            s._background_vector = \
             nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        return s._background_vector


    @property 
    def stats(self):
        try: 
            return self._stats
        except AttributeError: 
            pass
        self._stats = self.FeatureStats(self)
        return self._stats


    class FeatureStats:
        """Stores statistics about a L{FeatureScores}
        
        @ivar feats_total: Total number of features in the feature space.
        @ivar feats_masked: Number of masked features (unused, given zero score).
        @ivar feats_used: Number of features in use (non-zero scores).
        
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
            m = fs.mask
            
            s.feats_total = len(fs)
            s.feats_masked = sum(m)
            s.feats_used = s.feats_total - s.feats_masked
            
            s.pos_docs = fs.pdocs
            s.neg_docs = fs.ndocs
            
            s.pos_occurrences = int(nx.sum(fs.pos_counts[~m]))
            s.neg_occurrences = int(nx.sum(fs.neg_counts[~m]))
            
            s.pos_average = s.pos_occurrences / fs.pdocs if fs.pdocs > 0 else 0.0
            s.neg_average = s.neg_occurrences / fs.ndocs if fs.ndocs > 0 else 0.0
            
            s.pos_distinct = sum((fs.pos_counts != 0) & ~m)
            s.neg_distinct = sum((fs.neg_counts != 0) & ~m)


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
        # Mask to exclude zero-frequency features
        u = (docfreq_t != 0) 
        # Inverse Document Frequency (log N/df_t)
        idf = nx.log(N / docfreq_t[u])
        # Term frequency
        tf = (self.pos_counts[u] / nx.sum(self.pos_counts))
        # Calculate TF.IDF
        self._tfidf[u] = tf * idf
        return self._tfidf


    def get_best_tfidfs(self, count):
        """Construct a table about the features with the best TF.IDF.
        
        @param count: Maximum of rows to return.
        
        @return: List of (Feature ID, TFIDF, (feature, featureclass), 
        feature score, positive count, negative count)."""
        from heapq import nlargest
        best_tfidfs = nlargest(
            count, enumerate(self.tfidf), key=lambda x:x[1])
        return [ (t, tfidf, self.featmap[t], self.scores[t], 
                  self.pos_counts[t], self.neg_counts[t])
                  for t, tfidf in best_tfidfs ]


    def write_csv(self, stream, maxfeats=None):
        """Write features scores as CSV to an output stream.
        @param stream: Object supporting a write method.
        @param maxfeats: Optionally write only the top N features."""
        stream.write(u"score,positives,negatives,"\
                     u"numerator,denominator,"\
                     u"termid,type,term\n")
        s = self
        s.tfidf
        if maxfeats is None:
            maxfeats = len(s.scores)
        nmasked = sum(s.mask)
        logging.info("There are %d features, %d masked and %d unmasked",
                     len(s.mask), nmasked, len(s.mask)-nmasked)
        # Don't print more than maxfeats features (save space)
        numfeats = 0
        # Sort by decreasing score
        for (t, score) in sorted(
            enumerate(s.scores), key=lambda x:x[1], reverse=True):
            if s.mask[t]:
                continue
            numfeats += 1
            if numfeats > maxfeats:
                break
            stream.write(
                u'%.3f,%d,%d, %.2e,%.2e, %d,%s,"%s"\n' % 
                (s.scores[t], s.pos_counts[t], s.neg_counts[t], 
                 s.pfreqs[t], s.nfreqs[t], 
                 t, s.featmap[t][1], s.featmap[t][0]))



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