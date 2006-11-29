#!env python

from medline import *
from path import path
import tempfile
import unittest

class MedlineCacheTests(unittest.TestCase):
    """
    Tests MedlineCache updateCacheFromDir (so also makeDBEnv, putArticleList)
    """
    def setUp( self ):
        self.home = path(tempfile.mkdtemp(prefix="medline-"))

    def tearDown( self ):
        self.home.rmtree(ignore_errors=True)

    def test( self ):
        import test_xmlparse
        import xmlparse
        h = self.home
        xml = h/"test.xml"
        pmids = h/"pmids.txt"
        artdb = h/"articles.db"
        featdb = h/"features.db"
        fmap = FeatureMapping()
        m = MedlineCache(fmap,
                         xmlparse.ArticleParser(),
                         h,
                         artdb,
                         featdb,
                         h/"articles.txt",
                         h/"termcounts.pickle",
                         h/"processed.txt",
                         use_transactions=True,)
        xml.write_text(test_xmlparse.xmltext)
        m.updateCacheFromDir(h, save_delay=1)
        pmids.write_lines([ "1", "2" ])
        from article import getArticles
        a = getArticles(artdb, pmids)
        self.assertEqual(a[0].pmid, 1)
        self.assertEqual(a[1].pmid, 2)
        self.assertEqual(
            m.termcounts,
            {0: 2, 1: 2, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 1, 8: 1, 9: 1})
        self.assertEqual(
            fmap.feats,
            [(u'T6', 'mesh'), (u'T7', 'mesh'), (u'T4', 'mesh'),
             (u'T5', 'mesh'), (u'T2', 'mesh'), (u'T3', 'mesh'),
             (u'T1', 'mesh'), (u'0301-4851', 'issn'), (u'C2', 'mesh'),
             (u'C1', 'mesh')])

if __name__ == "__main__":
    unittest.main()
