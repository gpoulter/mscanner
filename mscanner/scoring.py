"""Calculate term and document scores

FeatureInfo -- Feature score calculation
iterScores() -- Iterate over document scores
filterDocuments() -- Filter iterScores() output for top N
partitionList() -- Split list into train/test partitions
retrievalTest() -- Count gold standard docs as a function of rank

                                   
"""

from __future__ import division

__license__ = """This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import logging as log
import numpy as nx

from mscanner import statusfile
from mscanner.support.utils  import selfupdate

class FeatureInfo:
    """Class to calculate and store all information about the
    features, their distribution, and their scores
    
    @note: FeatureInfo models a table with columns (pos_counts, neg_counts,
    pfreqs, nfreqs, mask, scores), and links to a FeatureMapping table (id,
    name, type, count) for Medline occurrence counts.

    Passed via constructor:
    
        @ivar pos_counts: Array of feature counts in positive documents. 
        
        @ivar neg_counts: Array of feature counts in negatives documents
        
        @ivar pdocs: Number of positive documents
        
        @ivar ndocs: Number of negative documents
    
        @ivar pseudocount: Either a float which is the pseudocount for all
        features, or numpy array of per-feature counts, or None which triggers
        once-off creation of per-features counts equal to Medline frequency if
        Bayesian calculation is in use.
        
        @ivar mask: Boolean array for masking some features scores to zero
        
        @ivar daniel: If true, use the JAMIA paper's smoothing heuristic
        of 10^-8 for terms found in positive but not negative (and visa
        versa).
        
        @ivar cutoff: If true, any features whose positive count is zero but
        nonetheless obtain a positive score (due to pseudocount priors), has
        its score set to zero.
    
    Added by constructor:
    
        @ivar getFrequencies: Method to calculate pfreq, nfreq
        
        @ivar postmask: Method for custom masking after scores are known
        
        @ivar featmap: featuremap.FeatureMapping object
    
    Set by getFeatureScores() and updateFeatureScores():
    
        @ivar pfreqs: Probability of each feature given positive article
        
        @ivar pfreqs: Probability of each feature given negative article
        
        @ivar scores: Array of feature scores

    """

    def __init__(self, 
                 featmap, 
                 pos_counts=None, 
                 neg_counts=None, 
                 pdocs=None, 
                 ndocs=None,
                 pseudocount=None, 
                 mask=None,
                 daniel=False, 
                 cutoff=False):
        """Initialise FeatureInfo object (parameters are instance variables)
        """
        if daniel:
            getFrequencies = self.getFrequenciesRubin05
        else:
            getFrequencies = self.getFrequenciesBayes
        if cutoff:
            postmask = self.maskRarePositives
        else:
            postmask = None
        selfupdate()
        
    def __len__(self):
        """Length is the number of features"""
        return len(self.featmap)
        
    def updateFeatureScores(self, pos_counts, neg_counts, pdocs, ndocs):
        """Change the feature counts and number of documents, clear
        old score calculations, and calculate new scores."""
        selfupdate()
        try:
            del self.scores
            del self.stats
        except AttributeError:
            pass
        return self.getFeatureScores()

    def getFeatureScores(self):
        """Calculate feature probabilities, and then calculate log likelihood
        support score of each feature.
        
        @return: Array of feature scores
        """
        pfreqs, nfreqs = self.getFrequencies()
        scores = nx.log(pfreqs / nfreqs)
        if self.mask is not None:
            pfreqs[self.mask] = 0
            nfreqs[self.mask] = 0
            scores[self.mask] = 0
        selfupdate()
        if self.postmask:
            self.postmask()
        return scores

    def getFrequenciesRubin05(s):
        """Uses Daniel Rubin's frequency smoothing heuristic, which
        replaces zero-frequency with 1e-8"""
        pfreqs = s.pos_counts / float(s.pdocs)
        nfreqs = s.neg_counts / float(s.ndocs)
        pfreqs[pfreqs == 0.0] = 1e-8
        nfreqs[nfreqs == 0.0] = 1e-8
        return pfreqs, nfreqs

    def getFrequenciesBayes(s):
        """Use Bayesian pseudocount vector, or a psuedocount float
         as a zero offset"""
        if s.pseudocount is None:
            s.pseudocount = nx.array(s.featmap.counts, nx.float32) / s.featmap.numdocs
        pfreqs = (s.pos_counts+s.pseudocount) / (2*s.pseudocount+s.pdocs)
        nfreqs = (s.neg_counts+s.pseudocount) / (2*s.pseudocount+s.ndocs)
        return pfreqs, nfreqs

    def maskRarePositives(s):
        """Removes positive-scoring features which fail to occur in the
        positive set"""
        s.scores[(s.scores > 0) & (s.pos_counts == 0)] = 0
        
    def maskLessThanThree(s):
        """Remove features occurring 2 or fewer times in total"""
        s.scores[(s.pos_counts + s.neg_counts) <= 2] = 0

    def getFeatureStats(self):
        """Recalculate statistics about the feature database
        
        Returns a Storage dictionary with the following keys:
    
        @ivar pos_occurrences: Total feature occurrences in positives
    
        @ivar neg_occurrences: Total feature occurrences in negatives
    
        @ivar feats_per_pos: Number of features per positive article
    
        @ivar feats_per_neg: Number of features per negative article
    
        @ivar distinct_feats: Number of distinct features
    
        @ivar pos_distinct_feats: Number of of distinct features in positives
    
        @ivar neg_distinct_feats: Number of of distinct features in negatives
        
        """
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
        return s

    def writeScoresBinary(self, stream):
        """Write feature scores as binary representation to stream"""
        stream.write(self.scores.tostring())

    def writeScoresCSV(self, stream):
        """Write features scores as CSV to an output stream"""
        if not hasattr(self, "scores"):
            self.getFeatureScores()
        stream.write(u"score,positives,negatives,pseudocount,"\
                     u"numerator,denominator,termid,type,term\n")
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
                u'%.3f,%d,%d,%.2e,%.2e,%.2e,%d,%s,"%s"\n' % 
                (s.scores[t], s.pos_counts[t], s.neg_counts[t], 
                 pseudocount[t], s.pfreqs[t], s.nfreqs[t],
                 t, s.featmap[t][1], s.featmap[t][0]))

def iterScores(docs, featscores, exclude=None):
    """Iterate over docs, yielding scores (calculated using featscores array),
    while skipping over members of exclude, and making status updates
    every 100,000 articles
    
    @param docs: Iterator over (integer doc ID, array of feature ID) pairs

    @param featscores: Array of feature scores (mapping feature ID to score)
    
    @param exclude: PMIDs to exclude from scoring
    
    @return: Iteration of (score, PMID) pairs
    """
    log.debug("Calculating article scores")
    marker = 0
    for idx, (docid, features) in enumerate(docs):
        if idx == marker:
            statusfile.update(idx)
            marker += 100000
        if exclude is None or docid not in exclude:
            yield nx.sum(featscores[features]), docid
    statusfile.update(None)    

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

def filterDocuments(scores, limit, threshold=None):
    """Return scores for documents given features and feature scores

    @param scores: Iterator over (score, pmid) pairs

    @param limit: Maximum number of results to return.

    @param threshold: Cutoff score for including an article in the results

    @return: List of (score, PMID) pairs in decreasing order of score
    """
    import heapq
    results = [(-100000, 0)] * limit
    ndocs = 0
    for score, docid in scores:
        if threshold is None or score >= threshold:
            ndocs += 1
            if score >= results[0][0]:
                heapq.heapreplace(results, (score,docid))
    if ndocs < limit:
        limit = ndocs
    return heapq.nlargest(limit, results)

def partitionList(pmids, subset_size=0.1, randomise=True):
    """Split list of objects into into two partitions
    
    @param pmids: List of PubMed IDs
    
    @param testsize: Amount to use for training (float for proportion,
    int for exact number)
    
    @param randomise: Set to False to disable random for debugging
    
    @return: Set of training IDs, and set of testing IDs
    """
    if randomise:
        import random
        random.shuffle(pmids)
    if isinstance(subset_size, float):
        subset_size = int(subset_size * len(pmids))
    return set(pmids[:subset_size]), set(pmids[subset_size:])

def retrievalTest(results, test):
    """Returns number of 

    @param results: List of result PubMed IDs
    
    @param test: Set of gold-standard articles to look for in query results

    @return: Array with of cumulative number of gold-standard articles at each 
    rank.
    """
    assert isinstance(test, set)
    TP_total = nx.zeros(len(results), nx.int32)
    for idx, pmid in enumerate(results):
        TP_total[idx] = TP_total[idx-1] if idx > 0 else 0
        if pmid in test:
            TP_total[idx] += 1
    return TP_total
