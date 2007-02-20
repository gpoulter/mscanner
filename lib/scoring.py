"""Calculate term and document scores, and print results

@author: Graham Poulter
                                   

FeatureScoreInfo -- Contains all information about feature scores
calculateFeatureScores() -- Calculate scores of features
filterDocuments() -- Return list of (docid,score) paris given documents and term scores
writeReport() -- Collate results of a query (calls all of above)
"""

from __future__ import division
from itertools import chain, izip
from path import path
import codecs
import heapq
import templates
import article
import numpy
import time

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
                 pseudocount, featmap, exclude_types=None, daniel=False):
        """Parameters are as for calculateFeatureScores, and also saved as instance variables.

        @param featmap: Mapping between Feature ID and (feature string, feature type)

        @param exclude_types: Types of features to exclude. Features of this
        type are or'd to the mask, where the default mask is those feature which
        have pos_counts[i] and neg_counts[i] equal to zero.
 
        """
        self.pos_counts = pos_counts
        self.neg_counts = neg_counts
        self.pdocs = pdocs
        self.ndocs = ndocs
        self.pseudocount = pseudocount
        self.featmap = featmap
        self.exclude_types = exclude_types
        self.daniel = daniel
        self.num_feats = len(self.featmap)
        self.mask = featmap.featureTypeMask(self.exclude_types)
        self.recalculateFeatureScores()
        self.recalculateFeatureStats()

    def recalculateFeatureScores(self):
        """Recalculate scores for features, and return the array of scores"""
        self.scores, self.pfreqs, self.nfreqs = calculateFeatureScores(
            self.pos_counts, self.neg_counts,
            self.pdocs, self.ndocs,
            self.pseudocount, self.mask, self.daniel)
        return self.scores

    def recalculateFeatureStats(self):
        """Recalculate statistics about the feature database"""
        self.pos_occurrences = int(numpy.sum(self.pos_counts)) 
        self.feats_per_pos = 0.0
        if self.pdocs > 0:
            self.feats_per_pos = self.pos_occurrences / self.pdocs 
        self.neg_occurrences = int(numpy.sum(self.neg_counts))
        self.feats_per_neg = 0.0
        if self.ndocs > 0:
            self.feats_per_neg = self.neg_occurrences / self.ndocs 
        self.pos_distinct_feats = len(numpy.nonzero(self.pos_counts)[0]) 
        self.neg_distinct_feats = len(numpy.nonzero(self.neg_counts)[0]) 

    def writeFeatureScoresCSV(self, stream):
        """Write features scores as CSV to an output stream"""
        stream.write(u"score,numerator,denominator,positives,negatives,termid,type,term\n")
        _ = self
        for line in sorted(
            [ (_.scores[t], _.pfreqs[t], _.nfreqs[t],
               _.pos_counts[t], _.neg_counts[t],
               t, _.featmap[t][1], _.featmap[t][0])
              for t in xrange(_.num_feats) if _.mask is None or not _.mask[t] ],
            key=lambda x:x[0], reverse=True):
            stream.write(u'%.3f,%.2e,%.2e,%d,%d,%d,%s,"%s"\n' % line)

def calculateFeatureScores(
    pos_counts, neg_counts,
    pdocs, ndocs,
    pseudocount, mask=None, daniel=False):
    """Return feature support scores based on relative frequency in positives vs negatives

    @param pos_counts: Array of feature counts in positive documents. 
    @param neg_counts: Array of feature counts in negatives documents

    @param pdocs: Number of positive documents
    @param ndocs: Number of negative documents

    @param pseudocount: Pseudocount for each term, e.g. if 0.1, we add 0.1
    occurrences across of the term, and add 0.2 (2*0.1) to the number of
    articles so that P(term|positive)+P(~term|positive)=1

    @param mask: Booolean array specifiying features to mask to zero

    @param daniel: If true, use the JAMIA paper's smoothing heuristic
    of 10^-8 for terms found in positive but not negative (and visa
    versa).

    @return: Array of feature scores (zeros for masked features),
    array of feature frequencies in positives, and array of
    feature frequencies of negatives.
    """
    num_feats = len(pos_counts) # number of distinct feature
    if daniel:
        pfreqs = pos_counts / float(pdocs)
        nfreqs = neg_counts / float(ndocs)
        pfreqs[pfreqs == 0.0] = 1e-8
        nfreqs[nfreqs == 0.0] = 1e-8
    else:
        pfreqs = (pos_counts+pseudocount) / (pdocs+2*pseudocount)
        nfreqs = (neg_counts+pseudocount) / (ndocs+2*pseudocount)
    scores = numpy.log(pfreqs / nfreqs)
    if mask is not None:
        scores[mask] = 0
    ## Uncomment to exclude "prior scores" on features we haven't encountered
    #dmask = (pos_counts + neg_counts == 0)
    #pfreqs[dmask] = 0
    #nfreqs[dmask] = 0
    #scores[dmask] = 0
    return scores, pfreqs, nfreqs

def filterDocuments(docs, featscores, limit=10000, threshold=0.0, statfile=lambda x:x):
    """Return scores for documents given features and feature scores

    @param docs: Iterable over (doc ID, array of feature ID) pairs

    @param featscores: Array of feature scores (mapping feature ID to score)

    @param limit: Maximum number of results to return.

    @param threshold: Cutoff score for including an article in the results

    @param statfile: Callable for posting status updates

    @return: Generator over (docid,score) pairs
    """
    results = [(-100000, 0)] * limit
    ndocs = 0
    marker = 0
    for idx, (docid, features) in enumerate(docs):
        if idx == marker:
            statfile(idx)
            marker += 100000
        score = numpy.sum(featscores[features])
        if score >= threshold:
            #print idx, docid, score
            ndocs += 1
            if score >= results[0][0]:
                heapq.heapreplace(results, (score,int(docid)))
    statfile(None)
    if ndocs < limit:
        limit = ndocs
    return [(pmid,score) for score,pmid in heapq.nlargest(limit, results)]

def writeReport(input, output, feature_info, configuration, artdb):
    """Write a report using the results of the classifier
    @param input: Iterator over (PubMed ID, Score)
    @param output: Iterator over (PubMed ID, Score)
    @param feature_info: A FeatureScoreInfo object with information about features
    @param configuration: Configuration module with remaining parameters
    @param artdb: Mapping from Pubmed ID to Article object
    """
    c = configuration
    rd = c.reportdir
    isinstance(feature_info, FeatureScoreInfo)
    feature_info.writeFeatureScoresCSV(codecs.open(rd/c.term_scores_name,"wb","utf-8"))
    # Citations for input articles
    templates.citations.run(
        dict(mode="input", 
             scores=sorted(input, key=lambda x:x[1], reverse=True), 
             articles=artdb),
        outputfile=file(rd/c.input_citations, "w"))
    # Citations for output articles
    templates.citations.run(
        dict(mode="output", 
             scores=sorted(output, key=lambda x:x[1], reverse=True), 
             articles=artdb),
        outputfile=file(rd/c.result_citations, "w"))
    # Index file
    templates.results.run(dict(
        time = time.strftime("%Y-%m-%d %H:%M:%S"),
        c = configuration,
        f = feature_info,
        num_results = len(output),
        lowest_score = output[-1][1],
        ), outputfile=(rd/c.index_file).open("w"))
    c.stylesheet.copy(rd/"style.css")
