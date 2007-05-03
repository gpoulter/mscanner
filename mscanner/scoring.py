"""Calculate term and document scores

calculateFeatureScores() -- Calculate scores of features
FeatureScoreInfo -- Contains all information about feature scores
filterDocuments() -- Return list of (docid,score) pairs given documents and term scores

                                   
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

import codecs
import heapq
import logging as log
import numpy as nx
import os
from path import path
import time

from mscanner import statusfile
from mscanner.featuremap import FeatureMapping
from mscanner.storage import Storage
from mscanner.gcheetah import FileTransaction, TemplateMapper

class FeatureScoreInfo:
    """Class to calculate and store all information about the
    features, their distribution, and their scores
    
    @ivar scores: Array of feature scores

    @ivar mask: Boolean array to mask out features of types in exclude_feats
    
    @ivar pos_occurrences: Total feature occurrences in positives

    @ivar neg_occurrences: Total feature occurrences in negatives

    @ivar feats_per_pos: Number of features per positive article

    @ivar feats_per_neg: Number of features per negative article

    @ivar distinct_feats: Number of distinct features

    @ivar pos_distinct_feats: Number of of distinct features in positives

    @ivar neg_distinct_feats: Number of of distinct features in negatives
    
    """

    def __init__(self, pos_counts, neg_counts, pdocs, ndocs,
                 pseudocount, featmap, exclude_types=None, 
                 daniel=False, cutoff=False):
        """Parameters are as for calculateFeatureScores, and are also kept as
        instance variables.

        @param featmap: Mapping between Feature ID and (feature string, feature
        type)

        @param exclude_types: Types of features to exclude (their scores are set
        to zero).
        
        @param pseudocount: Either a float which is the pseudocount for any
        feature. If "not pseudocount" (false/0/None), case we create an array of
        per-feature pseudocount, with each being the Medline frequency of the
        feature.
        
        @param pos_counts, neg_count, pdocs, ndocs, daniel, cutoff: Passed on
        direct to calculateFeatureScores.
  
        """
        isinstance(featmap, FeatureMapping)
        self.pos_counts = pos_counts
        self.neg_counts = neg_counts
        self.pdocs = pdocs
        self.ndocs = ndocs
        self.pseudocount = pseudocount
        if not pseudocount:
            self.pseudocount = nx.array(featmap.counts, nx.float32) / featmap.numdocs
        self.featmap = featmap
        self.exclude_types = exclude_types
        self.daniel = daniel
        self.cutoff = cutoff
        self.num_feats = len(self.featmap)
        self.mask = self.featmap.featureTypeMask(self.exclude_types)
        self.recalculateFeatureScores()
        self.recalculateFeatureStats()

    def recalculateFeatureScores(self):
        """Recalculate scores for features, and return the array of scores"""
        self.scores, self.pfreqs, self.nfreqs = calculateFeatureScores(
            self.pos_counts, self.neg_counts,
            self.pdocs, self.ndocs,
            self.pseudocount, self.mask, self.daniel, self.cutoff)
        return self.scores

    def recalculateFeatureStats(self):
        """Recalculate statistics about the feature database"""
        self.pos_occurrences = int(nx.sum(self.pos_counts)) 
        self.feats_per_pos = 0.0
        if self.pdocs > 0:
            self.feats_per_pos = self.pos_occurrences / self.pdocs 
        self.neg_occurrences = int(nx.sum(self.neg_counts))
        self.feats_per_neg = 0.0
        if self.ndocs > 0:
            self.feats_per_neg = self.neg_occurrences / self.ndocs 
        self.pos_distinct_feats = len(nx.nonzero(self.pos_counts)[0]) 
        self.neg_distinct_feats = len(nx.nonzero(self.neg_counts)[0]) 

    def writeFeatureScoresCSV(self, stream):
        """Write features scores as CSV to an output stream"""
        stream.write(u"score,positives,negatives,pseudocount,numerator,denominator,termid,type,term\n")
        _ = self
        if isinstance(self.pseudocount, float):
            pseudocount = nx.zeros_like(_.scores) + self.pseudocount
        else:
            pseudocount = self.pseudocount
        for t, score in sorted(enumerate(_.scores), key=lambda x:x[1], reverse=True):
            if _.mask and _.mask[t]:
                continue
            stream.write(
                u'%.3f,%d,%d,%.2e,%.2e,%.2e,%d,%s,"%s"\n' % 
                (_.scores[t], _.pos_counts[t], _.neg_counts[t], 
                 pseudocount[t], _.pfreqs[t], _.nfreqs[t],
                 t, _.featmap[t][1], _.featmap[t][0]))

def calculateFeatureScores(
    pos_counts, neg_counts, pdocs, ndocs, pseudocount, 
    mask=None, daniel=False, cutoff=False):
    """Return feature support scores based on relative frequency in positives vs
    negatives

    @param pos_counts: Array of feature counts in positive documents. 
    @param neg_counts: Array of feature counts in negatives documents

    @param pdocs: Number of positive documents
    @param ndocs: Number of negative documents

    @param pseudocount: Float with the global pseudocount for
    any term, e.g. if 0.1, we add 0.1 occurrences across of the term, and add
    0.2 (2*0.1) to the number of articles so that
    P(term|positive)+P(~term|positive)=1.  May instead be an array of floats
    specifying a different pseudocount for each feature.

    @param mask: Optional boolean array specifiying features to mask to zero

    @param daniel: If true, use the JAMIA paper's smoothing heuristic
    of 10^-8 for terms found in positive but not negative (and visa
    versa).
    
    @param cutoff: If true, any features whose positive count is zero but
    nonetheless obtain a positive score (due to pseudocount priors), has its
    score set to zero.
    
    @note: Features where pos_counts[i] and neg_counts[i] are both zero 
    have scores set to zero (disabled in favour of cutoff).

    @return: Array of feature scores (zeros for masked features),
    array of feature frequencies in positives, and array of
    feature frequencies of negatives.
    """
    # Daniel = Rubin2005 smoothing heuristic
    if daniel:
        pfreqs = pos_counts / float(pdocs)
        nfreqs = neg_counts / float(ndocs)
        pfreqs[pfreqs == 0.0] = 1e-8
        nfreqs[nfreqs == 0.0] = 1e-8
    # Use Bayesian pseudocount vector, or zero-offset psuedocount float
    else:
        pfreqs = (pos_counts+pseudocount) / (2*pseudocount+pdocs)
        nfreqs = (neg_counts+pseudocount) / (2*pseudocount+ndocs)
    # Log likelihood scores
    scores = nx.log(pfreqs / nfreqs)
    # Cutoff removes positive-scoring rare features
    if cutoff:
        scores[(scores > 0) & (pos_counts == 0)] = 0
    # Remove masked features from consideration
    if mask:
        pfreqs[mask] = 0
        nfreqs[mask] = 0
        scores[mask] = 0
    # Set "prior scores" for not-encountered features to zero?
    if False:
        dmask = (pos_counts + neg_counts == 0)
        pfreqs[dmask] = 0
        nfreqs[dmask] = 0
        scores[dmask] = 0
    return scores, pfreqs, nfreqs

def filterDocuments(docs, featscores, exclude=[], limit=10000, threshold=0.0):
    """Return scores for documents given features and feature scores

    @param docs: Iterator over (integer doc ID, array of feature ID) pairs

    @param featscores: Array of feature scores (mapping feature ID to score)
    
    @param exclude: Set of doc IDs to exclude from scoring

    @param limit: Maximum number of results to return.

    @param threshold: Cutoff score for including an article in the results

    @return: Generator over (docid,score) pairs
    """
    results = [(-100000, 0)] * limit
    ndocs = 0
    marker = 0
    for idx, (docid, features) in enumerate(docs):
        if docid in exclude:
            continue
        if idx == marker:
            statusfile.update(idx)
            marker += 100000
        score = nx.sum(featscores[features])
        if score >= threshold:
            #print idx, docid, score
            ndocs += 1
            if score >= results[0][0]:
                heapq.heapreplace(results, (score,docid))
    statusfile.update(None)
    if ndocs < limit:
        limit = ndocs
    return [(docid,score) for score,docid in heapq.nlargest(limit, results)]