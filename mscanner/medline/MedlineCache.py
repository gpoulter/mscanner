"""For updating the databases of articles and features"""

from __future__ import with_statement
from bsddb import db
import gzip
import logging as log
from path import path

from mscanner.medline import Shelf
from mscanner.medline.Article import Article
from mscanner.medline.FeatureDatabase import FeatureDatabase, FeatureStream
from mscanner.medline.FileTracker import FileTracker
from mscanner.medline.FeatureMapping import FeatureMapping


                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


class MedlineCache:    
    """Class for updating the Article DB, FeatureMapping, FeatureDatabase,
    FeatureStream, PMID list, and FileTracker.
    
    @ivar featmap: A FeatureMapping object for mapping string features to IDs
    @ivar db_env_home: Path to DB home directory 
    @ivar article_db: Path to article database
    @ivar feature_db: Path to feature database
    @ivar feature_stream: Path to feature stream file
    @ivar article_list: Path to list of article PMIDs
    @ivar narticles_path: Path to file containing the total number of PMIDs
    @ivar processed_path: Path to list of processed files
    @ivar use_transactions: If false, disable transaction engine
    """

    def __init__(
        self,
        featmap,
        db_env_home,
        article_db,
        feature_db,
        feature_stream,
        article_list,
        processed_path,
        narticles_path,
        use_transactions=True):
        """Constructor parameters set corresponding instance variables."""
        self.db_env_home = db_env_home
        self.featmap = featmap
        self.article_db = article_db
        self.feature_db = feature_db
        self.feature_stream = feature_stream
        self.article_list = article_list
        self.processed_path = processed_path
        self.narticles_path = narticles_path
        self.use_transactions = use_transactions
        self.recover = False


    def create_dbenv(self):
        """Create a Berkeley DB environment for transactions
        
        @return: DBEnv instance"""
        if not self.db_env_home.isdir():
            self.db_env_home.mkdir()
        dbenv = db.DBEnv()
        dbenv.set_lg_max(128*1024*1024) # 128Mb log files
        dbenv.set_tx_max(1) # 1 transaction at a time
        dbenv.set_cachesize(0, 8*1024*1024) # 8Mb shared cache
        flags = db.DB_INIT_MPOOL|db.DB_CREATE
        if self.use_transactions:
            flags |= db.DB_INIT_TXN
        if self.recover:
            flags |= db.DB_RECOVER # might use db.DB_RECOVER_FATAL
        dbenv.open(self.db_env_home, flags)
        return dbenv


    def add_articles(self, articles, dbenv):
        """Store Articles and feature lists in the databases
        
        @note: Databases are opened and closed inside each call, so that the
        user can Ctrl-C during the timed delay between files without corrupting
        the database.  Using transactions has too much overhead in time,
        and space used by the log files.
        
        @param articles: Iterator over Article objects
        
        @param dbenv: Database environment to use
        """
        log.info("Starting transaction to add articles")
        txn = None
        if self.use_transactions:
            txn = dbenv.txn_begin()
        try:
            artdb = Shelf.open(self.article_db, dbenv=dbenv, txn=txn)
            meshfeatdb = FeatureDatabase(self.feature_db, dbenv=dbenv, txn=txn)
            featstream = FeatureStream(open(self.feature_stream,"ab"))
            if not self.narticles_path.isfile():
                narticles = len(meshfeatdb)
            else:
                narticles = int(self.narticles_path.text())
            pmidlist = []
            for art in articles:
                # Refuse to add duplicates
                if art.pmid in meshfeatdb: 
                    continue
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
                featids = self.featmap.add_article(mesh=headings, qual=quals, issn=issns)
                meshfeatdb.setitem(art.pmid, featids, txn)
                featstream.write(art.pmid, featids)
            artdb.close()
            meshfeatdb.close()
            featstream.close()
            self.article_list.write_lines(pmidlist, append=True)
            self.narticles_path.write_text(str(narticles+len(pmidlist))+"\n")
            self.featmap.dump()
            if txn is not None:
                txn.commit()
        except Exception, e:
            if txn is not None:
                log.error("Aborting Transaction: Error %s", e)
                txn.abort()
            raise
        else:
            if txn is not None:
                log.info("Committed transaction")
            return len(pmidlist)


    def add_directory(self, medlinedir, save_delay=5):
        """Adds articles from XML files to MScanner's databases
        
        @param medlinedir: Path to a directory containing .xml.gz
        files
        
        @param save_delay: Pause this many seconds between calls to
        L{add_articles}"""
        import time
        filenames = medlinedir.files("*.xml") + medlinedir.files("*.xml.gz")
        tracker = FileTracker(self.processed_path)
        toprocess = tracker.toprocess(filenames)
        dbenv = self.create_dbenv()
        for idx, filename in enumerate(toprocess):
            log.info("Adding to cache: file %d out of %d (%s)", idx+1, len(toprocess), filename.name)
            for t in xrange(save_delay):
                log.debug("Saving in %d seconds...", save_delay-t)
                time.sleep(1)
            log.debug("Parsing XML file %s", filename.basename())
            try:
                if filename.endswith(".gz"):
                    infile = gzip.open(filename, 'r')
                else:
                    infile = open(filename, 'r')
                numadded = self.add_articles(Article.parse_medline_xml(infile), dbenv)
            finally:
                infile.close()
            log.debug("Added %d articles", numadded)
            tracker.add(filename)
            tracker.dump()
            log.info("Completed file %d out of %d (%s)", idx+1, len(toprocess), filename.name)
        dbenv.close()






