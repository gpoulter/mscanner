"""Manage Article and feature database

@author: Graham Poulter
                                    

"""

import os
from bsddb import db
import logging as log
from array import array
from path import path

import dbshelve
from featuredb import FeatureDatabase
from article import Article, FileTracker, FeatureMapping

class MedlineCache:
    """Manages the database of medline abstracts.  Function is to
    update the cache with new articles, and retrieve articles from the
    cache.
    """

    def __init__(
        self,
        featmap,
        parser,
        db_env_home,
        article_db,
        feature_db,
        article_list,
        processed_path,
        use_transactions=True):
        """Initialse a cache of the results of parsing medline.

        @param featmap: A FeatureMapping object for mapping string features to IDs
        @param parser: An ArticleParser object to parse XML
        @param db_env_home: Path to DB home directory 
        @param article_db: Path to article database
        @param feature_db: Path to feature database
        @param article_list: Path to list of article PMIDs
        @param processed_path: Path to list of processed files
        @param use_transactions: If false, disable transaction engine
        """
        self.db_env_home = db_env_home
        self.featmap = featmap
        self.parser = parser
        self.article_db = article_db
        self.feature_db = feature_db
        self.article_list = article_list
        self.processed_path = processed_path
        self.use_transactions = use_transactions
        self.recover = False

    def makeDBEnv(self):
        """Initialise DB environment for transactions"""
        if not self.db_env_home.isdir():
            self.db_env_home.mkdir()
        dbenv = db.DBEnv()
        dbenv.set_lg_max(512*1024*1024) # 512Mb log files
        dbenv.set_tx_max(1) # 1 transaction at a time
        dbenv.set_cachesize(0, 8*1024*1024) # 8Mb shared cache
        flags = db.DB_INIT_MPOOL|db.DB_CREATE
        if self.use_transactions:
            flags |= db.DB_INIT_TXN
        if self.recover:
            flags |= db.DB_RECOVER # db.DB_RECOVER_FATAL
        dbenv.open(self.db_env_home, flags)
        return dbenv

    def putArticleList(self, articles, dbenv):
        """Store Article objects and feature lists the databases

        @param articles: List of Article objects
        @param dbenv: Database environment to use
        """
        log.info("Starting transaction to add articles")
        txn = None
        if self.use_transactions:
            txn = dbenv.txn_begin()
        try:
            artdb = dbshelve.open(self.article_db, dbenv=dbenv, txn=txn)
            meshfeatdb = FeatureDatabase(self.feature_db, dbenv=dbenv, txn=txn)
            pmidlist = []
            for art in articles:
                # Refuse to add duplicates
                if art.pmid in meshfeatdb: continue
                # Store article, adding it to list of documents
                artdb[str(art.pmid)] = art
                pmidlist.append(str(art.pmid))
                # Get MeSH headings, qualifiers and ISSN from article
                headings = list()
                quals = list()
                issns = list()
                for term in art.meshterms:
                    headings.append(term[0])
                    if(len(term)>1):
                        for q in term[1:]:
                            if q not in quals:
                                quals.append(q)
                if art.issn:
                    issns = [art.issn]
                # Add features to feature mapping
                featids = self.featmap.addArticle(mesh=headings, qual=quals, issn=issns)
                meshfeatdb.setitem(art.pmid, featids, txn)
            artdb.close()
            meshfeatdb.close()
            self.article_list.write_lines(pmidlist, append=True)
            self.featmap.dump()
            if txn is not None:
                txn.commit()
        except Exception, e:
            raise
            if txn is not None:
                log.error("Aborting Transaction: Error %s", e)
                txn.abort()
        else:
            if txn is not None:
                log.info("Committed transaction")
            
    def updateCacheFromDir(self, medlinedir, save_delay=5):
        """Updates the cache given that medlinedir contains .xml.gz
        file to add to the cache and that we should save the inverse
        document after processing each savesteps files."""
        import time
        filenames = medlinedir.files("*.xml") + medlinedir.files("*.xml.gz")
        tracker = FileTracker(self.processed_path)
        toprocess = tracker.toprocess(filenames)
        dbenv = self.makeDBEnv()
        for idx, f in enumerate(toprocess):
            log.info("Adding to cache: file %d out of %d (%s)", idx+1, len(toprocess), f.name)
            for t in xrange(save_delay):
                log.debug("Saving in %d seconds...", save_delay-t)
                time.sleep(1)
            articles = list(self.parser.parseFile(f))
            log.debug("Parsed %d articles", len(articles))
            self.putArticleList(articles, dbenv)
            tracker.add(f)
            tracker.dump()
            log.info("Completed file %d out of %d (%s)", idx+1, len(toprocess), f.name)
        dbenv.close()
