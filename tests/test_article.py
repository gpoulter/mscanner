#!env python

import unittest
from path import path
from article import *

class FileTrackerTests(unittest.TestCase):
    """Test for FileTracker class

    Tests add, toprocess, dump, __init__
    """
    def setUp(self):
        self.home = path( '/tmp/test_filetracker' )
        self.home.rmtree( ignore_errors=True )
        self.home.mkdir()
    def tearDown( self ):
        self.home.rmtree( ignore_errors=True )
    def test( self ):
        t = FileTracker( self.home / 'test.txt' )
        t.add( path( "hack/a.xml" ) )
        t.add( path( "cough/b.xml" ) )
        self.assertEqual( t.toprocess([ path("foo/a.xml"), path("blah/c.xml") ]), ["blah/c.xml"] )
        t.dump()
        del t
        t = FileTracker( self.home / 'test.txt' )
        self.assertEqual( t, set( [ 'a.xml', 'b.xml' ] ) )
        
class FeatureMappingTests(unittest.TestCase):
    """Test for FeatureMapping class

    Tests: getids, getterms, __get__, dump, load
    """
    def test(self):
        fn = path( "/tmp/test_featuremapping" )
        fm = FeatureMapping( fn )
        self.assertEqual( fm.getids( ["A","B"] ), array("H",[0,1]) )
        self.assertEqual( fm.getterms( [0,1] ), [ "A", "B" ] )
        self.assertEqual( fm[0], "A" )
        self.assertEqual( fm[1], "B" )
        fm.dump()
        fm.load()
        fm.dump()
        fm.load()
        self.assertEqual( fm.term, ["A","B"] )
        self.assertEqual( fm.termid, {"A":0,"B":1} )
        try:
            fn.remove()
        except os.error:
            pass
        
class TermCountsTests(unittest.TestCase):
    """Test for TermCounts class

    Tests: add, __get__, dump, load, subtract
    """
    def test(self):
        t = TermCounts()
        t.add( [1,3] )
        t.add( [2,3] )
        self.assertEqual( t[1], 1 )
        self.assertEqual( t[2], 1 )
        self.assertEqual( t[3], 2 )
        self.assertEqual( t.docs, 2 )
        self.assertEqual( t.total, 4 )
        TermCounts.dump( t, '/tmp/test_termcounts' )
        t = TermCounts.load( '/tmp/test_termcounts' )
        s = TermCounts()
        s.add( [1,3] )
        r = t.subtract( s )
        self.assertEqual( r[1], 0 )
        self.assertEqual( r[2], 1 )
        self.assertEqual( r[3], 1 )
        self.assertEqual( r.docs, 1 )
        self.assertEqual( r.total, 2 )

class ArticleTests(unittest.TestCase):
    """Tests for article module functions

    Tests: chooseRandomLines, readPMIDFile, getArticles
    """
    def test(self):
        pass

if __name__ == "__main__":
    unittest.main()
