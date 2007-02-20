from path import path
import article
import scoring
import tempfile
import unittest
import numpy

class ScoringTests(unittest.TestCase):
    """Tests for scoring module functions"""
    def setUp(self):
        self.prefix = path(tempfile.mkdtemp(prefix="scoring-"))
        print self.prefix

    def tearDown(self):
        return
        self.prefix.rmtree(ignore_errors=True)
        
    def test_FeatureScoreInfo(self):
        """Implicitly tests calculateFeatureScores"""
        pfreqs = numpy.array([1,2,0])
        nfreqs = numpy.array([2,1,0])
        pdocs = 2
        ndocs = 3
        fm = article.FeatureMapping()
        f = scoring.FeatureScoreInfo(pfreqs, nfreqs, pdocs, ndocs, 0.1, fm)
        # With masking of unknown features
        #self.assert_(numpy.allclose(
        #    f.scores, numpy.array([-0.27193372,  1.02132061,  0.0 ])))
        # Without masking of unknown features
        self.assert_(numpy.allclose(
            f.scores, numpy.array([-0.27193372,  1.02132061,  0.37469345])))

    def test_filterDocuments(self):
        docs = { 1:[0,2], 2:[1,2], 3:[0,1] }
        fscores = numpy.array([1.0, 10.0, 100.0])
        self.assertEqual(numpy.sum(fscores[docs[1]]), 101.0)
        self.assertEqual(scoring.filterDocuments(docs.iteritems(), fscores, 2, 0 ), [ (2,110.0), (1,101.0) ])
        self.assertEqual(scoring.filterDocuments(docs.iteritems(), fscores, 10, 102.0 ), [ (2,110.0) ])

    def test_writeReport(self):
        pfreqs = numpy.array([1,1,2])
        nfreqs = numpy.array([1,1,0])
        pdocs = 2
        ndocs = 1
        articles = {
            "1111": article.Article(1111,"T","A",meshterms=set(["A","B"])),
            "2222": article.Article(2222,"T","A",meshterms=set(["A","C"])),
            "3333": article.Article(3333,"T","A",meshterms=set(["B","C"])),
            }
        scores = [(3333,3.5), (2222,2.0), (1111,1.0)]
        featmap = [("A","T"), ("B","T"), ("C","Q")]

if __name__ == "__main__":
    unittest.main()
