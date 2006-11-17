#!env python

from medline import *
import os
from path import path
import unittest
import warnings
warnings.filterwarnings("ignore")

class FeatureDatabaseTests(unittest.TestCase):
    """
    Tests FeatureDatabase getitem, __iter__, iteritems, in, len
    """
    def test( self ):
        d = FeatureDatabase(value_type="H")
        d.setitem( 1, array("H",[1,3]) )
        d.setitem( 2, array("H",[2,3]) )
        self.assertEqual( d.getitem(1), array("H",[1,3]) )
        self.assertEqual( d.getitem(2), array("H",[2,3]) )
        self.assertRaises( KeyError, d.getitem, 3 )
        self.failUnless( 1 in d )
        self.failUnless( 2 in d )
        self.failIf( 3 in d )
        self.assertEqual( d.keys(), [2,1] )
        self.assertEqual( list( d.__iter__() ), [2,1] )
        self.assertEqual( list( d.iteritems() ), [ (2,array("H",[2,3])), (1,array("H",[1,3])) ] )
        self.assertEqual( len(d), 2 )
        d.delitem(2)
        self.failIf( 2 in d )

class MedlineCacheTests(unittest.TestCase):
    """
    Tests MedlineCache updateCacheFromDir (so also makeDBEnv, putArticleList)
    """
    def setUp( self ):
        self.home = path(os.tempnam())
        self.home.mkdir()

    def tearDown( self ):
        self.home.rmtree(ignore_errors=True)

    def test( self ):
        import test_xmlparse
        import xmlparse
        h = self.home
        m = MedlineCache( FeatureMapping(),
                          xmlparse.ArticleParser(),
                          h,
                          h/"articles.db",
                          h/"features.db",
                          h/"articles.txt",
                          h/"termcounts.pickle",
                          h/"processed.txt" )
        (h/"test.xml").write_text( test_xmlparse.xmltext )
        m.updateCacheFromDir( h, save_delay=1 )
        (h/"pmids.xml").write_lines( [ "1", "2" ] )
        from article import getArticles
        a = getArticles( h/"articles.db", h/"pmids.xml" )
        self.assertEqual( a[0].pmid, 1 )
        self.assertEqual( a[1].pmid, 2 )
        self.assertEqual( m.termcounts, {0: 2, 1: 2, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2} )

if __name__ == "__main__":
    unittest.main()
