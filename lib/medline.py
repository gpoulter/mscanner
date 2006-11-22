"""Manage Article and feature database

@author: Graham Poulter
                                    

FeatureDatabase -- Persistent mapping from document ID to array of Feature IDs
MedlineCache -- Parses and adds articles to the databases

"""

import os
import struct
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
        @param flags: Opening flags (r,rw,w,c,n)
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
        """Return an array object of values for a given key"""
        values_packed = self.db.get(struct.pack("i",int(key)), txn=txn)
        if values_packed is None:
            raise KeyError("Record %d not found in feature database" % key)
        return array(self.value_type, values_packed)

    def __getitem__(self, key):
        return self.getitem(key)

    def setitem(self, key, values, txn=None ):
        """Associate integer key with an array object of values"""
        if not isinstance( values, array ):
            raise TypeError("values must be an array('%s')" % self.value_type)
        self.db.put(struct.pack("i",int(key)), values.tostring(), txn=txn)
        
    def delitem(self, key, txn=None ):
        """Delete a given key from the database"""
        self.db.delete(struct.pack("i",int(key)), txn=txn)

    def __len__(self):
        return len(self.db)

    def __contains__(self, key):
        return self.db.has_key(struct.pack("i",int(key)))

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
        article_list_path,
        termcounts_path,
        processed_path ):
        """Initialse a cache of the results of parsing medline.

        @param meshdb: A FeatureMapping object for MeSH terms <-> 16-bit IDs
        @param parser: An ArticleParser object to parse XML
        @param db_env_home: Path to DB home directory 
        @param article_db_path: Path to article database
        @param feature_db_path: Path to feature database
        @param article_list_path: Path to list of article PMIDs
        @param termcounts_path: Path to the TermCounts pickle
        @param processed_path: Path to list of processed files 
        """
        self.db_env_home = db_env_home
        self.meshdb = meshdb
        self.parser = parser
        self.article_db_path = article_db_path
        self.feature_db_path = feature_db_path
        self.article_list_path = article_list_path
        self.termcounts_path = termcounts_path
        self.processed_path = processed_path
        self.termcounts = TermCounts.load( termcounts_path )

    def makeDBEnv( self ):
        """Initialise DB environment for transactions"""
        try: self.db_env_home.mkdir()
        except os.error: pass
        dbenv = db.DBEnv()
        dbenv.set_lg_max(512*1024*1024) # 512Mb log files
        dbenv.set_tx_max(1) # 1 transaction at a time
        dbenv.set_cachesize(0, 8*1024*1024) # 8Mb shared cache
        dbenv.open(self.db_env_home, db.DB_INIT_MPOOL|db.DB_INIT_TXN|db.DB_CREATE)
        return dbenv

    def putArticleList(self, articles, dbenv):
        """Write a list of Article objects to the cache"""
        # Starting transaction
        log.info("Starting transaction to add articles")
        txn = dbenv.txn_begin()
        try:
            artdb = dbshelve.open(self.article_db_path, dbenv=dbenv, txn=txn)
            featdb = FeatureDatabase(self.feature_db_path, dbenv=dbenv, txn=txn)
            artlist = file( self.article_list_path, "a" )
            for art in articles:
                # Refuse to add or overwrite duplicates
                if not art.pmid in featdb:
                    termids = self.meshdb.getids(art.meshterms)
                    artdb[str(art.pmid)] = art
                    featdb.setitem(art.pmid, termids, txn)
                    artlist.write("%d\n" % art.pmid)
                    self.termcounts.add(termids)
            artdb.close()
            featdb.close()
            artlist.close()
            txn.commit()
            self.meshdb.dump()
            TermCounts.dump(self.termcounts, self.termcounts_path)
        except Exception, e:
            log.error( "Aborting Transaction: Error %s", e )
            txn.abort()
            raise
        else:
            log.info( "Successfully committed transaction to add articles" )
            
    def updateCacheFromDir(self, medlinedir, save_delay=5):
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
