#!/usr/bin/env python

"""Manage Article and feature database

@author: Graham Poulter
                                    

FeatureDatabase -- Persistent mapping from document ID to array of Feature IDs
MedlineCache -- Parses and adds articles to the databases

"""

import os
import struct
import unittest
import logging as log
from bsddb import db
from array import array
from path import path

import dbshelve
from article import Article, FileTracker, FeatureMapping, TermCounts

class FeatureDatabase:
    """Persistent mapping from integer key to array objects of numerical values"""
    
    def __init__(self, filename=None, flags='c', mode=0660, dbenv=None, txn=None, value_type='H'):
        """Initialise database

        @param filename: Database file
        @param flags: Opening flats (r,rw,w,c,n)
        @param mode: Numeric file permissions
        @param dbenv: Optional database environment
        @param txn: Optional database transaction
        @param value_type: Typecode for packing/unpacking of value structs, typically 'H' or 'i'
        """
        self.db = db.DB( dbenv )
        self.value_type = value_type
        if isinstance(flags, basestring):
            if flags == 'r':
                flags = db.DB_RDONLY
            elif flags == 'rw':
                flags = 0
            elif flags == 'w':
                flags = db.DB_CREATE
            elif flags == 'c':
                flags = db.DB_CREATE
            elif flags == 'n':
                flags = db.DB_TRUNCATE | db.DB_CREATE
            else:
                raise db.DBError("Flag %s is not in 'r', 'rw', 'w', 'c' or 'n'"  % str(flags))
        self.db.open( filename, None, db.DB_HASH, flags, mode, txn=txn )

    def __del__(self):
        if hasattr(self,"db"):
            self.close()

    def close(self):
        """Close the database.  Do not use this object after doing so"""
        self.db.close()
        delattr(self, "db")

    def getitem(self, key, txn=None):
        """Return array of values for a document id"""
        values_packed = self.db.get(struct.pack("i",key), txn=txn)
        if values_packed is None:
            raise KeyError("Record %d not found in feature database" % key)
        return array(self.value_type, values_packed)

    def __getitem__(self, key):
        return self.getitem(key)

    def setitem(self, key, values, txn=None ):
        """Set array of values for a given key"""
        if not isinstance( values, array ):
            raise TypeError("values must be an array('%s')" % self.value_type)
        self.db.put(struct.pack("i",key), values.tostring(), txn=txn)
        
    def delitem(self, key, txn=None ):
        """Delete a given key from the database"""
        self.db.delete(struct.pack("i",key), txn=txn)

    def __len__(self):
        return len(self.db)

    def __contains__(self, key):
        return self.db.has_key(struct.pack("i",key))

    def keys(self):
        return [ k for k in self ]

    def __iter__( self ):
        cur = self.db.cursor()
        rec = cur.first(dlen=0, doff=0)
        while rec is not None:
            yield struct.unpack("i",rec[0])[0]
            rec = cur.next(dlen=0, doff=0)
        cur.close()

    def iteritems( self ):
        cur = self.db.cursor()
        rec = cur.first()
        while rec is not None:
            yield struct.unpack("i",rec[0])[0], array(self.value_type,rec[1])
            rec = cur.next()
        cur.close()

class MedlineCache:
    """Manages the database of medline abstracts.  Function is to
    update the cache with new articles, and retrieve articles from the
    cache.

    @note: Other instance variables are parameters to __init__
    @ivar termcounts: Mapping from term ID to number of occurrences
    """

    def __init__(
        self,
        meshdb,
        parser,
        db_env_home,
        article_db_path,
        feature_db_path,
        termcounts_path,
        processed_path ):
        """Initialse a cache of the results of parsing medline.

        @param meshdb: A FeatureMapping object for MeSH terms <-> 16-bit IDs
        @param parser: An ArticleParser object to parse XML
        @param db_env_home: Path to DB home directory 
        @param article_db_path: Path to article database
        @param feature_db_path: Path to feature database
        @param termcounts_path: Path to the TermCounts pickle
        @param processed_path: Path to list of processed files 
        """
        self.db_env_home = db_env_home
        self.meshdb = meshdb
        self.parser = parser
        self.article_db_path = article_db_path
        self.feature_db_path = feature_db_path
        self.termcounts_path = termcounts_path
        self.processed_path = processed_path
        self.termcounts = TermCounts.load( termcounts_path )

    def makeDBEnv( self ):
        """Initialise DB environment for transactions"""
        try: self.db_env_home.mkdir()
        except os.error: pass
        dbenv = db.DBEnv()
        dbenv.set_lg_max( 512*1024*1024 ) # 512Mb log files
        dbenv.set_tx_max( 1 ) # 1 transaction at a time
        dbenv.set_cachesize( 0, 8*1024*1024 ) # 8Mb shared cache
        dbenv.open( self.db_env_home, db.DB_INIT_MPOOL|db.DB_INIT_TXN|db.DB_CREATE|db.DB_RECOVER_FATAL )
        return dbenv

    def putArticleList( self, articles, dbenv ):
        """Write a list of Article objects to the cache"""
        # Starting transaction
        log.info("Starting transaction to add articles")
        txn = dbenv.txn_begin()
        try:
            artdb = dbshelve.open( self.article_db_path, dbenv=dbenv, txn=txn )
            featdb = FeatureDatabase( self.feature_db_path, dbenv=dbenv, txn=txn )
            for art in articles:
                # Refuse to add or overwrite duplicates
                if not art.pmid in featdb:
                    termids = self.meshdb.getids( art.meshterms )
                    artdb[str(art.pmid)] = art
                    featdb.setitem( art.pmid, termids, txn )
                    self.termcounts.add( termids )
            artdb.close()
            featdb.close()
            txn.commit()
            self.meshdb.dump()
            TermCounts.dump( self.termcounts, self.termcounts_path )
        except Exception, e:
            log.error( "Aborting Transaction: Error %s", e )
            txn.abort()
            raise
        else:
            log.info( "Successfully committed transaction to add articles" )
            
    def updateCacheFromDir( self, medlinedir, save_delay=5 ):
        """Updates the cache given that medlinedir contains .xml.gz
        file to add to the cache and that we should save the inverse
        document after processing each savesteps files."""
        import time
        filenames = medlinedir.files( "*.xml" ) + medlinedir.files( "*.xml.gz" )
        tracker = FileTracker( self.processed_path )
        toprocess = tracker.toprocess( filenames )
        dbenv = self.makeDBEnv()
        for idx, f in enumerate( toprocess ):
            log.info( "Adding to cache: file %d out of %d (%s)", idx+1, len(toprocess), f.name )
            for t in xrange(save_delay):
                log.debug( "Saving in %d seconds...", save_delay-t )
                time.sleep(1)
            articles = list( self.parser.parseFile( f ) )
            self.putArticleList( articles, dbenv )
            tracker.add(f)
            tracker.dump()
            log.info( "Completed file %d out of %d (%s)", idx+1, len(toprocess), f.name )
        dbenv.close()

class _FeatureDatabaseTests(unittest.TestCase):
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

class _MedlineCacheTests(unittest.TestCase):
    def setUp( self ):
        self.home = h = path( '/tmp/medline_test' )
        h.rmtree( ignore_errors=True )
        try: h.mkdir()
        except os.error: pass
    def tearDown( self ):
        self.home.rmtree( ignore_errors=True )
    def test( self ):
        import xmlparse
        h = self.home
        m = MedlineCache( FeatureMapping(),
                          xmlparse.ArticleParser(),
                          h,
                          h/"articles.db",
                          h/"features.db",
                          h/"termcounts.pickle",
                          h/"processed.txt" )
        (h/"test.xml").write_text( xmlparse.xmltext )
        m.updateCacheFromDir( h, save_delay=1 )
        (h/"pmids.xml").write_lines( [ "1", "2" ] )
        from article import getArticles
        a = getArticles( h/"articles.db", h/"pmids.xml" )
        self.assertEqual( a[0].pmid, 1 )
        self.assertEqual( a[1].pmid, 2 )
        self.assertEqual( m.termcounts, {0: 2, 1: 2, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2} )

if __name__ == "__main__":
    unittest.main()
