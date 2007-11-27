"""Synchronised updating of databases with incoming Medline records."""

from __future__ import with_statement
from bsddb import db
import gzip
import logging
from path import path
import time

from mscanner.configuration import rc
from mscanner.medline.Article import Article
from mscanner.medline.Databases import ArticleData, FeatureData
from mscanner.medline.FeatureStream import Date2Integer


                                     
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



class Updater:    
    """Updates the databases with incoming Medline records. Reads data from XML
    files, and writes to L{ArticleData}, L{FeatureData} and L{FileTracker}.
    
    @ivar adata: L{ArticleData} to update article databases.
    @ivar fdata_mesh: L{FeatureData} to update MeSH article features.
    @ivar fdata_all: L{FeatureData} to update all article features.
    @ivar tracker: Path to list of processed files.
    @ivar stopwords: Set of words not to use as features.
    """


    def __init__(self, adata, fdata_mesh, fdata_all, tracker, stopwords):
        """Constructor parameters set corresponding instance variables.
        @param stopwords: Path to file of stopwords.
        @param tracker: Path to file to use with L{FileTracker}.
        """
        self.adata = adata
        self.fdata_mesh = fdata_mesh
        self.fdata_all = fdata_all
        self.tracker = tracker
        self.stopwords = [] if stopwords is None else stopwords.lines(retain=False)


    @staticmethod
    def Defaults():
        """Construct an Updater using the RC parameter defaults, and 
        a database environment."""
        dbenv = Updater.dbenv(rc.db_env_home)
        return Updater(ArticleData.Defaults(dbenv), 
                       FeatureData.Defaults_MeSH(dbenv),
                       FeatureData.Defaults_All(dbenv), 
                       rc.tracker, rc.stopwords)


    @staticmethod
    def dbenv(envdir, recover=False):
        """Create a Berkeley DB environment (for transactions mostly).
        @param envdir: Path to DB environment directory.
        @return: DBEnv instance"""
        if not envdir.isdir():
            envdir.mkdir()
        dbenv = db.DBEnv()
        dbenv.set_lg_max(128*1024*1024) # 128Mb log files
        dbenv.set_tx_max(1) # 1 transaction at a time
        dbenv.set_cachesize(0, 8*1024*1024) # 8Mb shared cache
        flags = db.DB_INIT_MPOOL|db.DB_CREATE
        # Might try db.DB_RECOVER_FATAL
        if recover: flags |= db.DB_RECOVER 
        # Not using transactions due to time/space overhead.
        if False: flags |= db.DB_INIT_TXN
        dbenv.open(envdir, flags)
        return dbenv
    
    
    def add_articles(self, articles):
        """Update feature and article databases with a list of new article
        objects. Dispatches to the L{ArticleData} and L{FeatureData} instances.
        @param articles: List of Article objects """
        self.adata.add_articles(articles)
        self.fdata_mesh.add_articles((str(a.pmid), 
            Date2Integer(a.date_completed), 
            a.mesh_features()) for a in articles)
        self.fdata_all.add_articles((str(a.pmid), 
            Date2Integer(a.date_completed), 
            a.all_features(self.stopwords)) for a in articles)


    def add_directory(self, medline=None, save_delay=5):
        """Adds articles from XML files to MScanner's databases

        @param medline: Directory with .xml.gz files. If None we default to
        C{rc.medline}.
        
        @param save_delay: Seconds to pause between files.
        """
        if medline is None: medline = rc.medline
        # List of input files
        infiles = medline.files("*.xml") + medline.files("*.xml.gz")
        # Set of completed files
        if self.tracker.isfile():
            done = set(self.tracker.lines(retain=False))
        else:
            done = set()
        # Ordered list of non-completed files
        todo = sorted([f for f in infiles if f.basename() not in done])
        # Loop over non-completed files
        for idx, filename in enumerate(todo):
            logging.info("Adding to cache: file %d out of %d (%s)", 
                         idx+1, len(todo), filename.name)
            for t in xrange(save_delay):
                logging.debug("Saving in %d seconds...", save_delay-t)
                time.sleep(1)
            logging.debug("Parsing XML file %s", filename.basename())
            try:
                if filename.endswith(".gz"):
                    infile = gzip.open(filename, 'r')
                else:
                    infile = open(filename, 'r')
                articles = list(Article.parse_medline_xml(infile))
                self.add_articles(articles)
            finally:
                infile.close()
            # Update the file tracker
            done.add(filename.basename())
            self.tracker.write_lines(sorted(done)) # Don't break here!
            logging.info("Added %d articles from file %d out of %d (%s)", 
                         len(articles), idx+1, len(todo), filename.name)
            del articles

