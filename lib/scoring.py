#!/usr/bin/env python

"""Calculate term and document scores, and print results

@author: Graham Poulter
                                   

getTermScores() -- Calculate termid:termscore mapping, given term counts in pos/neg corpus
scoreDocument() -- Calculate score for document given term IDs and term scores
filterDocuments() -- Return list of (docid,score) paris given documents and term scores
writeTermScoresCSV() -- Output term scores to CSV file
writeTermScoresHTML() -- Output term scores to HTML file
writeResultScoresTXT() -- Output results (docid,score) to TXT file
writeCitationsHTML() -- Output citations for results to HTML file
writeReport() -- Collate results of a query (calls all of above)

"""

from __future__ import division
import codecs
import unittest
from path import path
import templates

def getTermScores( positive, negative, pseudocount=0.1, daniel=False ):
    """Return log-likelihood support scores for MeSH terms.

    @param positive: TermCounts instance for positive articles. 

    @param negative: TermCounts instance for negative articles. 

    @param pseudocount: Pseudocount for each term, e.g. if 0.1, we
    pretend that each term has an additional 0.1 occurrences across
    the articles, adding 2*0.1 to the number of articles so that
    P(term|positive)+P(~term|positive)=1

    @param daniel: If true, use the JAMIA paper's smoothing heuristic
    of 10^-8 for terms found in positive but not negative (and visa
    versa), which greatly exaggerates scores for such terms.

    @rtype: C{{termid:(float,float,float,int,int)}}

    @return: Mapping of term IDs to (support score, frequency in
    positive, negative corpus, and count in positive, negative corpus)
    """
    import math
    score = dict()
    def makeScore(termid,pcount,ncount):
        if daniel:
            if pcount > 0:
                pfreq = pcount / positive.docs
            else:
                pfreq = 1e-8
            if ncount > 0:
                nfreq = ncount / negative.docs
            else:
                nfreq = 1e-8
        else:
            pfreq = (pcount+pseudocount) / (positive.docs+2*pseudocount)
            nfreq = (ncount+pseudocount) / (negative.docs+2*pseudocount)
        score[termid] = (math.log(pfreq/nfreq),pfreq,nfreq,pcount,ncount)
    for (termid,pcount) in positive.iteritems():
        if termid in score: continue
        ncount = negative.get( termid, 0 )
        makeScore(termid,pcount,ncount)
    for (termid,ncount) in negative.iteritems():
        if termid in score: continue
        pcount = positive.get( termid, 0 )
        makeScore(termid,pcount,ncount)
    return score

def scoreDocument( features, feature_scores ):
    """Return document score given a feature list and feature scores
    @param features: List of feature IDs
    @param feature_scores: Mapping of feature ID to score
    """
    score = 0.0
    for f in features:
        if f in feature_scores: 
            score += feature_scores[f][0]
    return score

def filterDocuments( docs, feature_scores, limit=-100000.0 ):
    """Return scores for documents given features and feature scores
    @param docs: Iteratable over (doc ID, feature ID list) pairs
    @param feature_scores: Mapping from feature ID to score
    @param limit: If int, it is maximum number of results to return.
    If float, it is the cutoff score.
    @return: List of (score,docid) pairs
    """
    results = list()
    if isinstance( limit, float ):
        for docid, features in docs:
            score = scoreDocument(features,feature_scores)
            if score >= limit:
                results.append( (score,docid) )
        results.sort( reverse=True )
    elif isinstance( limit, int ):
        ndocs = 0
        results = [ (-100000,0) for i in range( limit ) ]
        from heapq import heapreplace
        for docid, features in docs:
            ndocs += 1
            score = scoreDocument(features,feature_scores)
            if score >= results[0][0]:
                heapreplace( results, (score,docid) )
        results.sort( reverse=True )
        if ndocs < limit: del results[ndocs:]
    return results

def writeTermScoresCSV(f, meshdb, scores, pfreqs, nfreqs):
    """Write term scores to CSV file"""
    LINE = "%.3f,%.3f,%.3f,%d,%d,%d,%s\n"
    f.write("score,numerator,denominator,positives,negatives,termid,term\n")
    for termid, s in sorted( scores.iteritems(), key=lambda x:x[1][0], reverse=True ):
        f.write(LINE % (s[0],s[1],s[2],s[3],s[4],termid,meshdb[termid]))

def writeTermScoresHTML(f, meshdb, scores, pfreqs, nfreqs, pseudocount):
    """Write term scores to HTML file"""
    templates.termscore.run(dict(
        numterms = len(meshdb),
        pseudocount = pseudocount,
        pdocs = pfreqs.docs,
        ndocs = nfreqs.docs,
        poccurs = pfreqs.total,
        noccurs = nfreqs.total,
        pterms = len(pfreqs),
        nterms = len(nfreqs),
        scores = scores,
        meshdb = meshdb
        ),
        outputfile=f)

def writeReport(
    scores,
    meshdb,
    termscores,
    pfreq,
    nfreq,
    prefix,
    stylesheet,
    pseudocount,
    limit,
    posfile,
    articles=None,
    ):
    """Write a report using the results of the classifier
    @param scores: List of (score,docid) pairs, in decreasing order of score
    @param meshdb: termid:term mapping (FeatureMapping instance)
    @param termscores: termid:(score,...) mapping
    @param pfreq: termid:count mappings for positives 
    @param nfreq: termid:count mappings for negatives
    @param prefix: Directory in which to place output
    @param stylesheet: Path to CSS stylesheet
    @param pseudocount: From @{getTermScores}
    @param limit: From @{filterDocuments}
    @param posfile: Path to positive docids
    @param articles: Mapping from str(pmid) to Article object
    """
    terms_csv = prefix/"termscores.csv"
    terms_html = prefix/"termscores.html"
    res_txt = prefix/"resultscores.txt"
    input_html = prefix/"input_citations.html"
    result_html= prefix/"result_citations.html"
    index_html = prefix/"index.html"
    # Term scores
    writeTermScoresCSV(file(terms_csv,'w'), meshdb, termscores, pfreq, nfreq)
    writeTermScoresHTML(file(terms_html,'w'), meshdb, termscores, pfreq, nfreq, pseudocount)
    # Result scores
    f = file(res_txt,'w')
    for score, pmid in scores:
        f.write( "%-13d%.5f\n" % ( pmid, score ) )
    # Citations for input articles
    import article
    posids = list(article.readPMIDFile(posfile))
    templates.citations.run(
        dict(mode="input", scores=[(0,p) for p in posids], articles=articles),
        outputfile=file(input_html, "w"))
    # Citations for output articles
    templates.citations.run(
        dict(mode="output", scores=scores, articles=articles),
        outputfile=file(result_html, "w"))
    # Index file
    templates.results.run(dict(
        pseudocount = pseudocount,
        limit = limit,
        num_results = len(scores),
        lowest_score = scores[-1][0],
        posfile = posfile.basename(),
        terms_csv = terms_csv.basename(),
        terms_html = terms_html.basename(),
        res_txt = res_txt.basename(),
        input_html = input_html.basename(),
        result_html = result_html.basename(),
        ),
        outputfile=index_html.open("w"))
    posfile.copy(prefix / posfile.basename())
    stylesheet.copy(prefix / "style.css")

class _ScoringTests(unittest.TestCase):

    class Freqs(dict):
        def __init__(self, docs, total, contents):
            dict.__init__(self)
            self.docs = docs
            self.total = total
            self.update( contents )
        
    def test_getTermScores(self):
        pfreqs = self.Freqs( 2, 3, { 100:1, 200:2 } )
        nfreqs = self.Freqs( 2, 3, { 100:2, 200:1 } )
        termscores = getTermScores( pfreqs, nfreqs, pseudocount=0 )
        self.assertEqual( termscores,
                          {100: (-0.69314718055994529, 0.5, 1.0, 1, 2),
                           200: (0.69314718055994529, 1.0, 0.5, 2, 1)}
                          )

    def test_filterDocuments(self):
        docs = { 1:["A","C"], 2:["B","C","D"] }
        feat_scores = { "A":(1,), "B":(10,), "C":(100,) }
        self.assertEqual( scoreDocument( docs[1], feat_scores ), 101.0 )
        self.assertEqual( scoreDocument( docs[2], feat_scores ), 110.0 )
        self.assertEqual( filterDocuments( docs.iteritems(), feat_scores, 1 ), [ (110.0,2) ] )
        self.assertEqual( filterDocuments( docs.iteritems(), feat_scores, 102.0 ), [ (110.0,2) ] )

    def test_writeReport(self):
        pfreqs = self.Freqs( 2, 3, { 1:1, 2:2 } )
        nfreqs = self.Freqs( 2, 3, { 1:2, 2:1 } )
        from article import Article
        articles = {
            "1111": Article(1111,"T","A",set(["A","B"])),
            "2222": Article(2222,"T","A",set(["A","C"])),
            "3333": Article(3333,"T","A",set(["B","C"])),
            }
        path("/tmp/positives.txt").write_text("2222\n3333\n")
        termscores = getTermScores( pfreqs, nfreqs, pseudocount=0 )
        writeReport(
            scores = [ (1.0,1111), (2.0,2222), (3.5,3333) ],
            meshdb = { 1:"A", 2:"B", 3:"C" },
            termscores = termscores,
            pfreq = pfreqs,
            nfreq = nfreqs,
            prefix = path("/tmp/result_"),
            stylesheet = path("templates/style.css"),
            pseudocount = 0,
            limit = 3,
            posfile = path("/tmp/positives.txt"),
            articles = articles,
            )

if __name__ == "__main__":
    unittest.main()
