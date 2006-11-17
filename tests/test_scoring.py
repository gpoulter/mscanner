#!env python

import os
from path import path
from scoring import *
import unittest
import warnings
warnings.filterwarnings("ignore")

class ScoringTests(unittest.TestCase):
    """Tests for scoring module functions

    Tests: getTermScores, filterDocuments, writeReport
    Implicit: scoreDocument, writeTermScoresCSV, writeTermScoresHTML
    """
    class Freqs(dict):
        def __init__(self, docs, total, contents):
            dict.__init__(self)
            self.docs = docs
            self.total = total
            self.update(contents)

    def setUp(self):
        self.prefix = path(os.tempnam())
        self.prefix.mkdir()

    def tearDown(self):
        self.prefix.rmtree(ignore_errors=True)
        
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
        self.assertEqual(scoreDocument( docs[1], feat_scores ), 101.0)
        self.assertEqual(scoreDocument( docs[2], feat_scores ), 110.0)
        self.assertEqual(filterDocuments( docs.iteritems(), feat_scores, 1, 0 ), [ (110.0,2) ])
        self.assertEqual(filterDocuments( docs.iteritems(), feat_scores, 10, 102.0 ), [ (110.0,2) ])

    def test_writeReport(self):
        pfreqs = self.Freqs( 2, 3, { 1:1, 2:2 } )
        nfreqs = self.Freqs( 2, 3, { 1:2, 2:1 } )
        from article import Article
        articles = {
            "1111": Article(1111,"T","A",set(["A","B"])),
            "2222": Article(2222,"T","A",set(["A","C"])),
            "3333": Article(3333,"T","A",set(["B","C"])),
            }
        (self.prefix/"positives.txt").write_text("2222\n3333\n")
        termscores = getTermScores( pfreqs, nfreqs, pseudocount=0 )
        writeReport(
            scores = [(1.0,1111), (2.0,2222), (3.5,3333)],
            meshdb = { 1:"A", 2:"B", 3:"C" },
            termscores = termscores,
            pfreq = pfreqs,
            nfreq = nfreqs,
            prefix = self.prefix,
            stylesheet = path("../lib/templates/style.css"),
            pseudocount = 0,
            limit = 3,
            posfile = self.prefix/"positives.txt",
            articles = articles,
            )

if __name__ == "__main__":
    unittest.main()
