"""Calculate term and document scores, and print results

@author: Graham Poulter
                                   

FeatureScores -- Class containing all information about feature scores
calculateFeatureScores() -- Calculate scores of features
filterDocuments() -- Return list of (docid,score) paris given documents and term scores
writeReport() -- Collate results of a query (calls all of above)
"""

from __future__ import division
from itertools import chain, izip
from path import path
import codecs
import templates
import article
import numpy

class FeatureScoreInfo:
    """Class to calculate and store all information about the
    features, their distribution, and their scores"""

    def __init__(self, pos_counts, neg_counts, pdocs, ndocs,
                 pseudocount, featmap, exclude_types=None, daniel=False):
        """Parameters are as for calculateFeatureScores, and also saved as instance variables.

        @param featmap: Mapping between Feature ID and (feature string, feature type)

        @param exclude_types: Types of features to exclude (mask is calculated from this and featmap)

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
        self.pos_counts = pos_counts
        self.neg_counts = neg_counts
        self.pdocs = pdocs
        self.ndocs = ndocs
        self.pseudocount = pseudocount
        self.featmap = featmap
        self.exclude_types = exclude_types
        self.daniel = daniel
        # Calculate masked features
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

    @param pseudocount: Pseudocount for each term, e.g. if 0.1, we
    pretend that each term has an additional 0.1 occurrences across
    the articles, adding 2*0.1 to the number of articles so that
    P(term|positive)+P(~term|positive)=1

    @param mask: Booolean array specifiying features to mask to zero

    @param daniel: If true, use the JAMIA paper's smoothing heuristic
    of 10^-8 for terms found in positive but not negative (and visa
    versa), which greatly exaggerates scores for such terms.

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
    return scores, pfreqs, nfreqs

def filterDocuments(docs, featscores, limit=10000, threshold=0.0, statfile=None):
    """Return scores for documents given features and feature scores

    @param docs: Iteratable over (doc ID, feature ID list) pairs

    @param featscores: Mapping from feature ID to score

    @param limit: Maximum number of results to return.

    @param threshold: Cutoff score for including an article in the results

    @param statfile: Name of file to post status updates

    @return: Generator over (docid,score) pairs
    """
    results = list()
    ndocs = 0
    results = [ (-100000, 0) for i in range(limit) ]
    from heapq import heapreplace
    for idx, (docid, features) in enumerate(docs):
        if statfile is not None and idx % 100000 == 0:
            statfile.update(idx)
        score = numpy.sum(featscores[features])
        if score >= threshold:
            #print idx, docid, score
            ndocs += 1
            if score >= results[0][0]:
                heapreplace(results, (score,int(docid)))
    if statfile is not None:
        statfile.update(None)
    results.sort(reverse=True)
    if ndocs < limit:
        del results[ndocs:]
    return numpy.array([pmid for score,pmid in results],dtype=numpy.int32), \
           numpy.array([score for score,pmid in results],dtype=numpy.float32)

def writeReport(
    input_pmids,
    input_scores,
    result_pmids,
    result_scores,
    feature_info,
    configuration,
    artdb):
    """Write a report using the results of the classifier
    @param input_pmids: Input PubMed IDs
    @param input_pmids: Input scores
    @param result_pmids: Result PubMed IDs
    @param result_scores: Result scores
    @param feature_info: A FeatureScoreInfo object with information about features
    @param configuration: Configuration module with remaining parameters
    @param artdb: Mapping from Pubmed ID to Article object
    """
    c = configuration
    rd = c.reportdir
    feature_info.writeFeatureScoresCSV(codecs.open(rd/c.term_scores_name,"wb","utf-8"))
    # Citations for input articles
    input_pairs = sorted(zip(input_pmids, input_scores), key=lambda x:x[1], reverse=True)
    templates.citations.run(
        dict(mode="input", scores=input_pairs, articles=artdb),
        outputfile=file(rd/c.input_citations, "w"))
    # Citations for output articles
    templates.citations.run(
        dict(mode="output", scores=zip(result_pmids, result_scores), articles=artdb),
        outputfile=file(rd/c.result_citations, "w"))
    # Index file
    templates.results.run(dict(
        time = time.strftime("%Y-%m-%d %H:%M:%S"),
        c = configuration,
        f = feature_info,
        num_results = len(result_scores),
        lowest_score = result_scores[-1],
        ), outputfile=(rd/c.index_file).open("w"))
    c.stylesheet.copy(rd/"style.css")
