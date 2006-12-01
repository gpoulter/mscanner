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
        fmap = FeatureMapping(h/"featuremap.txt")
        m = MedlineCache(fmap,
                         xmlparse.ArticleParser(),
                         h,
                         artdb,
                         featdb,
                         h/"articles.txt",
                         h/"processed.txt",
                         use_transactions=True,)
        xml.write_text(test_xmlparse.xmltext)
        m.updateCacheFromDir(h, save_delay=1)
        pmids.write_lines(["1", "2"])
        from article import getArticles
        a = getArticles(artdb, pmids)
        self.assertEqual(a[0].pmid, 1)
        self.assertEqual(a[1].pmid, 2)
        self.assertEqual(fmap.freqs, [3, 3, 3, 3, 3, 3, 3, 2])
        self.assertEqual(
            fmap.feats,
            [(u'T1', 'mesh'), (u'T2', 'mesh'), (u'T3', 'mesh'),
             (u'T6', 'mesh'), (u'Q4', 'qual'), (u'Q5', 'qual'),
             (u'Q7', 'qual'), (u'0301-4851', 'issn')])
            
if __name__ == "__main__":
    unittest.main()
