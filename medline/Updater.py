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
from mscanner.medline.FeatureStream import DateAsInteger


                                     
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
    files, and writes to L{ArticleData} and L{FeatureData} objects.
    
    @ivar adata: L{ArticleData} instance containing article database.
    @ivar fdata_list: List of L{FeatureData} instances for article representations.
    @ivar tracker: Path to the list of completed Medline XML files.
    @ivar stopwords: Set of words not to use as features.
    """


    def __init__(self, adata, fdata_list, tracker):
        """Constructor parameters set corresponding instance variables."""
        self.adata = adata
        self.fdata_list = fdata_list
        self.tracker = tracker


    def close(self):
        """Close all L{ArticleData} and L{FeatureData} instances"""
        self.adata.close()
        for fdata in self.fdata_list:
            fdata.close()


    def sync(self):
        """Synchronise the FeatureData objects to disk, which can take up to 10
        minutes if there are a lot of features in the FeatureMap."""
        logging.info("PLEASE WAIT while overwriting feature maps.")
        for fdata in self.fdata_list:
            fdata.sync()
        logging.info("DONE overwriting feature maps.")


    @staticmethod
    def Defaults(featurespaces):
        """Construct an Updater using default path names..
        
        @param featurespaces: List of (name,ftype) pairs, where name specifies a
        feature-getting method on L{Article}, and ftype is uint16/uint32 for the
        FeatureDatabase numpy type."""
        return Updater(
            ArticleData.Defaults(), 
            [FeatureData.Defaults(name,ftype,rdonly=False) 
             for name,ftype in featurespaces], 
            rc.articles_home/rc.tracker)
    

    def regenerate(self):
        """Regenerate L{FeatureData} from the L{ArticleData}."""
        self.adata.regenerate_artlist()
        for fdata in self.fdata_list:
            fdata.regenerate(self.adata.artdb.itervalues())


    def add_articles(self, articles, sync=True):
        """Add a list of new article objects to L{ArticleData} and
        L{FeatureData} instances.

        @param articles: List of Article objects to add.
        
        @param sync: If True, synchronise feature maps after adding."""
        logging.warn("Adding articles to databases. DO NOT INTERRUPT!")
        self.adata.add_articles(articles)
        for fdata in self.fdata_list:
            fdata.add_articles(articles)
        if sync: self.sync()


    def add_directory(self, medline, save_delay=5):
        """Adds articles from XML files to MScanner's databases.  
        
        @note: Do not CTRL-C the updater on Windows! We write the feature map
        at the end inside a finally. However, windows CTRL-C exits immediately
        (on Unix we get to clean up). There is no workaround - just don't break
        mid-way on windows.
        
        @param medline: Directory with .xml.gz files. If None we default to
        C{rc.medline}.
        
        @param save_delay: Seconds to pause between files.
        """
        # List of input files
        infiles = medline.files("*.xml") + medline.files("*.xml.gz")
        # Set of completed files
        if self.tracker.isfile():
            done = set(self.tracker.lines(retain=False))
        else:
            done = set()
        # Ordered list of non-completed files
        todo = sorted([f for f in infiles if f.basename() not in done])
        try:
            # Loop over non-completed files
            for idx, filename in enumerate(todo):
                for t in xrange(save_delay):
                    logging.debug("Pausing for %d seconds...", save_delay-t)
                    time.sleep(1)
                try:
                    logging.info("Parsing XML file %d out of %d (%s)", 
                                 idx+1, len(todo), filename.name)
                    if filename.endswith(".gz"):
                        infile = gzip.open(filename, 'r')
                    else:
                        infile = open(filename, 'r')
                    articles = list(Article.parse_medline_xml(infile))
                except KeyboardInterrupt:
                    logging.info("Safely interrupted.")
                    raise
                finally:
                    infile.close()
                try:
                    self.add_articles(articles, sync=False)
                    done.add(filename.basename())
                    self.tracker.write_lines(sorted(done))
                    logging.info("Added %d articles from file %d out of %d (%s)", 
                                 len(articles), idx+1, len(todo), filename.name)
                    del articles
                except KeyboardInterrupt:
                    logging.error("Unsafely interrupted!!!")
                    raise
        finally:
            self.sync()
