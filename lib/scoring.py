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
from path import path
from codecs import open
import templates
import article
import numpy

def getTermScores(
    positive,
    negative,
    pdocs,
    ndocs,
    pseudocount=0.01,
    daniel=False,
    featmap=None,
    exclude=None):
    """Return feature support scores based on relative frequency in positives vs negatives

    @param positive: Frequencies of features in positive documents. 
    @param negative: Frequencies of features in negatives documents

    @param pdocs: Number of positive documents
    @param ndocs: Number of negative documents

    @param pseudocount: Pseudocount for each term, e.g. if 0.1, we
    pretend that each term has an additional 0.1 occurrences across
    the articles, adding 2*0.1 to the number of articles so that
    P(term|positive)+P(~term|positive)=1

    @param daniel: If true, use the JAMIA paper's smoothing heuristic
    of 10^-8 for terms found in positive but not negative (and visa
    versa), which greatly exaggerates scores for such terms.

    @param featmap: Mapping from feature ID to (feature string, feature type)

    @param exclude: List of feature types to exclude from scoring

    @rtype: C{{termid:(float,float,float,int,int)}}

    @return: Mapping of term IDs to (support score, frequency in
    positive, negative corpus, and count in positive, negative corpus)
    """
    import math
    def makeScore(termid,pcount,ncount):
        if daniel:
            pfreq, nfreq = 1e-8, 1e-8
            if pcount > 0:
                pfreq = pcount / pdocs
            if ncount > 0:
                nfreq = ncount / ndocs
        else:
            pfreq = (pcount+pseudocount) / (pdocs+2*pseudocount)
            nfreq = (ncount+pseudocount) / (ndocs+2*pseudocount)
        feature_type = featmap[termid][1]
        if exclude is None or feature_type not in exclude:
            return [math.log(pfreq/nfreq), pfreq, nfreq, pcount, ncount]
        else:
            return [0, 0, 0, 0, 0]
    return numpy.array([makeScore(termid,pcount,negative[termid]) for termid, pcount in enumerate(positive)])

def scoreDocument(features, featscores):
    """Return document score given a feature list and feature scores
    @param features: Array of feature IDs
    @param featscores: Array of (score,pfreq,nfreq,pcount,ncount) for feature IDs
    """
    return numpy.sum(featscores[:,0][features])

def filterDocuments(docs, featscores, limit=10000, threshold=0.0, statfile=None):
    """Return scores for documents given features and feature scores

    @param docs: Iteratable over (doc ID, feature ID list) pairs

    @param featscores: Mapping from feature ID to score

    @param limit: Maximum number of results to return.

    @param threshold: Cutoff score for including an article in the results

    @param statfile: Name of file to post status updates

    @return: List of (score,docid) pairs
    """
    results = list()
    ndocs = 0
    results = [ (-100000, 0) for i in range(limit) ]
    from heapq import heapreplace
    for idx, (docid, features) in enumerate(docs):
        if statfile is not None and idx % 100000 == 0:
            statfile.update(idx)
        score = scoreDocument(features, featscores)
        if score >= threshold:
            ndocs += 1
            if score >= results[0][0]:
                heapreplace(results, (score,docid))
    if statfile is not None:
        statfile.update(None)
    results.sort(reverse=True)
    if ndocs < limit:
        del results[ndocs:]
    return results

def writeTermScoresCSV(f, featmap, featscores):
    """Write term scores to CSV file"""
    f.write(u"score,numerator,denominator,positives,negatives,termid,term,type\n")
    for s in sorted(
        [ (s[0],s[1],s[2],s[3],s[4],termid,featmap[termid][0],featmap[termid][1])
          for termid, s in enumerate(featscores) if s[3] != 0 or s[4] != 0 ],
        key=lambda x:x[0], reverse=True):
        f.write(u'%.3f,%.2e,%.2e,%d,%d,%d,"%s",%s\n' % s)

def writeTermScoresHTML(f, featmap, featscores, pdocs, ndocs, pseudocount):
    """Write term scores to HTML file"""
    pfreq, nfreq = 0.0, 0.0
    pfreqs, nfreqs = featscores[:,3], featscores[:,4]
    if pdocs:
        pfreq = numpy.sum(pfreqs)/pdocs
    if ndocs:
        nfreq = numpy.sum(nfreqs)/ndocs
    scorelines = sorted(
        [(s[0],s[1],s[2],s[3],s[4],termid,featmap[termid][0],featmap[termid][1])
         for termid, s in enumerate(featscores)],# if s[3] != 0 or s[4] != 0],
        key=lambda x:x[0], reverse=True)
    templates.termscore.run(dict(
        numterms = len(featmap),
        pseudocount = pseudocount,
        pdocs = pdocs,
        ndocs = ndocs,
        poccurs = int(numpy.sum(pfreqs)),
        noccurs = int(numpy.sum(nfreqs)),
        pterms = len(numpy.nonzero(pfreqs)[0]),
        nterms = len(numpy.nonzero(nfreqs)[0]),
        pfreq = pfreq,
        nfreq = nfreq,
        scores = scorelines
       ), outputfile=f)

def writeReport(
    scores,
    featmap,
    pdocs,
    ndocs,
    termscores,
    prefix,
    stylesheet,
    pseudocount,
    limit,
    posfile,
    articles=None,
   ):
    """Write a report using the results of the classifier
    @param scores: List of (score,docid) pairs, in decreasing order of score
    @param featmap: termid:term mapping (FeatureMapping instance)
    @param pdocs: Number of positive documents
    @param ndocs: Number of negative documents
    @param termscores: termid:(score,numerator,denominator,pfreq,nfreq) mapping
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
    writeTermScoresCSV(open(terms_csv,'w','utf-8'), featmap, termscores)
    writeTermScoresHTML(open(terms_html,'w','ascii','xmlcharrefreplace'),
                        featmap, termscores, pdocs, ndocs, pseudocount)
    # Result scores
    f = file(res_txt,'w')
    for score, pmid in scores:
        f.write("%-13s%.5f\n" % (pmid, score))
    # Citations for input articles
    import article
    posids = list(article.readPMIDFile(posfile, articles))
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
    if not (prefix / posfile.basename()).exists():
        posfile.copy(prefix / posfile.basename())
    stylesheet.copy(prefix / "style.css")
