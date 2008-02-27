"""Synchronised updating of databases with incoming Medline records."""

from __future__ import with_statement
from bsddb import db
import gzip
import logging
from path import path
import time

from mscanner.configuration import rc
from mscanner.medline.Article import Article
from mscanner.medline.FeatureData import FeatureData
from mscanner.medline.FeatureStream import DateAsInteger
from mscanner.medline import Shelf


                                     
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
    files, and writes to article database and L{FeatureData} indexes.
    
    @ivar artdb: L{Shelf} for looking up Article objects by PubMed ID
    @ivar fdata_list: List of L{FeatureData} instances for article representations.
    @ivar tracker: Path to the list of completed Medline XML files.
    @ivar stopwords: Set of words not to use as features.
    """


    def __init__(self, artdb, fdata_list, tracker):
        """Constructor parameters set corresponding instance variables."""
        self.artdb = artdb
        self.fdata_list = fdata_list
        self.tracker = tracker


    def close(self):
        """Close article database and L{FeatureData} instances"""
        self.artdb.close()
        for fdata in self.fdata_list:
            fdata.close()


    @staticmethod
    def Defaults(featurespaces):
        """Construct an Updater using default path names.
        @param featurespaces: List of feature space names to load, each of
        which specifies a feature-getting method on L{Article}."""
        # Create new database root if necessary
        base = rc.articles_home
        if not base.exists(): base.makedirs()
        return Updater(
            Shelf.open(base/rc.articledb),
            [FeatureData.Defaults(fs, rdonly=False) for fs in featurespaces], 
            rc.articles_home/rc.tracker)
    

    def regenerate(self):
        """Regenerate L{FeatureData}, possibly from the article database."""
        for fdata in self.fdata_list:
            fdata.regenerate(self.artdb)


    def add_articles(self, articles):
        """Add a list of new article objects to article database and L{FeatureData}.
        @param articles: List of Article objects to add.
        """
        logging.warn("Adding articles to databases. DO NOT INTERRUPT!")
        for art in articles:
            self.artdb[str(art.pmid)] = art
        self.artdb.sync()
        for fdata in self.fdata_list:
            fdata.add_articles(articles)


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
        # Loop over non-completed files
        for idx, filename in enumerate(todo):
            for t in xrange(save_delay):
                logging.debug("Pausing for %d seconds...", save_delay-t)
                time.sleep(1)
            # Parse the XML file
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
            # Add the articles to the databases
            try:
                self.add_articles(articles)
                done.add(filename.basename())
                self.tracker.write_lines(sorted(done))
                logging.info("Added %d articles from file %d out of %d (%s)", 
                             len(articles), idx+1, len(todo), filename.name)
                del articles
            except KeyboardInterrupt:
                logging.error("Unsafely interrupted!!!")
                raise
