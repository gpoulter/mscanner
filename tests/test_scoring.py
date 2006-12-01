#!env python

from path import path
from scoring import *
import tempfile
import unittest
import numpy

class ScoringTests(unittest.TestCase):
    """Tests for scoring module functions

    Tests: getTermScores, filterDocuments, writeReport
    Implicit: scoreDocument, writeTermScoresCSV, writeTermScoresHTML
    """
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="scoring-"))

    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)
        
    def test_getTermScores(self):
        pfreqs = numpy.array([1,2])
        nfreqs = numpy.array([2,1])
        pdocs = 2
        ndocs = 2
        termscores = getTermScores(pfreqs, nfreqs, pdocs, ndocs, pseudocount=0)
        self.assertEqual(termscores.tolist(), 
                          [[-0.69314718055994529, 0.5, 1.0, 1, 2],
                           [0.69314718055994529, 1.0, 0.5, 2, 1]])

    def test_filterDocuments(self):
        docs = { 1:[0,2], 2:[1,2] }
        feat_scores = numpy.array([(1,0,0,0,0),(10,0,0,0,0),(100,0,0,0,0)])
        self.assertEqual(scoreDocument(docs[1], feat_scores), 101.0)
        self.assertEqual(scoreDocument(docs[2], feat_scores), 110.0)
        self.assertEqual(filterDocuments(docs.iteritems(), feat_scores, 1, 0 ), [ (110.0,2) ])
        self.assertEqual(filterDocuments(docs.iteritems(), feat_scores, 10, 102.0 ), [ (110.0,2) ])

    def test_writeReport(self):
        pfreqs = numpy.array([1,1,2])
        nfreqs = numpy.array([1,1,0])
        pdocs = 2
        ndocs = 1
        from article import Article
        articles = {
            "1111": Article(1111,"T","A",set(["A","B"])),
            "2222": Article(2222,"T","A",set(["A","C"])),
            "3333": Article(3333,"T","A",set(["B","C"])),
            }
        (self.prefix/"positives.txt").write_text("2222\n3333\n")
        termscores = getTermScores(pfreqs, nfreqs, pdocs, ndocs, pseudocount=0.1)
        writeReport(
            scores = [(1.0,1111), (2.0,2222), (3.5,3333)],
            featmap = [ ("A","T"), ("B","T"), ("C","Q") ],
            pdocs = pdocs,
            ndocs = ndocs,
            termscores = termscores,
            prefix = self.prefix,
            stylesheet = path("../lib/templates/style.css"),
            pseudocount = 0.1,
            limit = 3,
            posfile = self.prefix/"positives.txt",
            articles = articles,
            )

if __name__ == "__main__":
    unittest.main()
