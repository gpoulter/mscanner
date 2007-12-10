"""Syncrhonising around the article and feature databases. They allow the
database to be loaded at the same time, and for articles are added
simultaneously to all of them databases."""

from __future__ import with_statement
from contextlib import closing
import logging
import numpy as nx

from mscanner.configuration import rc
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.medline.FeatureDatabase import FeatureDatabase
from mscanner.medline.FeatureStream import FeatureStream, Date2Integer
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


class FeatureData:
    """Wraps the files handling feature ID information - namely the
    L{FeatureMapping}, L{FeatureDatabase} and L{FeatureStream}.

    @ivar featmap: L{FeatureMapping} between feature names and feature IDs.
    
    @ivar featuredb: Mapping from PubMed ID to list of features.
    
    @ivar fstream: L{FeatureStream} of (PMID, date, feature vector).
    
    @ivar rdonly: Boolean for opening databases read-only.
    """
    
    def __init__(self, featmap, featdb, fstream, ftype, dbenv=None, rdonly=True):
        """Constructor. 
        @param dbenv: Optional Berkeley database environment"""
        self.rdonly = rdonly
        logging.debug("Loading feature mapping from %s", featmap.basename())
        self.featmap = FeatureMapping(featmap, ftype=ftype)
        self.featuredb = FeatureDatabase(featdb, "r" if rdonly else "c", 
                                         dbenv=dbenv, ftype=ftype)
        self.fstream = FeatureStream(fstream, ftype, rdonly)


    @staticmethod
    def Defaults_MeSH(dbenv=None, rdonly=True):
        """Construct using the RC defaults for MeSH feature space."""
        return FeatureData(rc.featuremap_mesh, rc.featuredb_mesh, 
                           rc.featurestream_mesh, nx.uint16, dbenv, rdonly)


    @staticmethod
    def Defaults_All(dbenv=None, rdonly=True):
        """Construct using RC defaults for the MeSH+Abstract feature space."""
        return FeatureData(rc.featuremap_all, rc.featuredb_all, 
                           rc.featurestream_all, nx.uint32, dbenv, rdonly)


    def close(self):
        """Shut down the databases"""
        if not self.rdonly:
            self.featmap.dump()
        self.featuredb.close()
        self.fstream.close()
        
        
    def sync(self):
        """Flush databases to disk"""
        if not self.rdonly:
            self.featmap.dump()
            self.featuredb.sync()
            self.fstream.flush()


    def add_articles(self, articles):
        """Incrementally add new articles to the existing feature 
        database, stream and feature map.  
        
        @param articles: Iterable of (PMID, date, features), as (string,
        YYYYMMDD integer, dictionary). The features dict maps string feature
        type to list of feature strings."""
        if self.rdonly:
            raise NotImplementedError("Attempt to write to read-only databases")
        for pmid, date, features in articles:
            if pmid not in self.featuredb:
                self.featmap.add_article(features)
                featvec = self.featmap.get_vector(features)
                self.featuredb[str(pmid)] = featvec
                self.fstream.additem(pmid, date, featvec)
        
        
    def regenerate(self, articles):
        """Regenerate feature map, feature stream and feature database. We only
        regenerate those files that are currently empty (because they were
        moved/deleted).
        
        @param articles: Iterable of (PMID, date, features), as (string,
        YYYYMMDD integer, dictionary). The features dict maps string feature
        type to list of feature strings."""
        if self.rdonly:
            raise NotImplementedError("Attempt to write to read-only databases")
        do_featmap = len(self.featmap) == 0
        do_stream = self.fstream.filename.size == 0
        do_featuredb = len(self.featuredb) == 0
        if not (do_featmap or do_stream or do_featuredb):
            logging.info("Not regenerating feature dbs as they are non-empty.")
            return
        if do_featmap: 
            logging.info("Regen feature map %s.", self.featmap.filename.basename())
        if do_stream: 
            logging.info("Regen feature stream %s.", self.fstream.filename.basename())
        if do_featuredb: 
            logging.info("Regen feature database %s.", self.featuredb.filename.basename())
        for count, (pmid, date, features) in enumerate(articles):
            if count % 10000 == 0:
                logging.debug("Regenerated %d articles.", count)
            if do_featmap:
                self.featmap.add_article(features)
            if do_stream or do_featuredb:
                featvec = self.featmap.get_vector(features)
                if do_stream:
                    self.fstream.additem(pmid, date, featvec)                
                if do_featuredb:
                    self.featuredb[str(pmid)] = featvec
        self.sync()



class ArticleData:
    """Wraps the files containing Article data - namely the Article database,
    and the list of PubMed IDs and dates.  
    
    @ivar artdb: Mapping from PubMed ID to Article object.
    
    @ivar artlist_path: Path to file listing PubMed ID and record date.
    """
    
    def __init__(self, artdb, artlist, dbenv=None):
        """Constructor for setting attributes to be used by the remaining
        methods."""
        self.artdb = Shelf.open(artdb, dbenv=dbenv)
        self.artlist_path = artlist


    @staticmethod
    def Defaults(dbenv=None):
        """Construct an instance using default RC parameters."""
        return ArticleData(rc.articledb, rc.articlelist, dbenv)


    def close(self):
        """Closes the article databases"""
        self.artdb.close()


    @property
    def article_list(self):
        """Array of PubMed IDs in the database. We read these from
        L{artlist_path}.  Takes a few minutes to load the first time!"""
        try:
            return self._article_list
        except AttributeError: 
            logging.info("Loading article list property.")
            if self.artlist_path.exists():
                with open(self.artlist_path, "rb") as f:
                    self._article_list = nx.array([int(t.split()[0]) for t in f], nx.uint32)
            else:
                self._article_list = nx.array([], nx.uint32)
            return self._article_list


    def regenerate_artlist(self):
        """Regenerate the list of articles from the article database."""
        if len(self.article_list) != 0:
            logging.info("Not regenerating article list as it is non-empty.")
            return
        del self._article_list
        with open(self.artlist_path, "wb") as f:
            for idx, (pmid, art) in enumerate(self.artdb.iteritems()):
                if idx % 10000 == 0:
                    logging.debug("Wrote %d to the article list." % idx)
                f.write("%s %08d\n" % (pmid, Date2Integer(art.date_completed)))


    def add_articles(self, articles):
        """Add new articles to both the database and the article list.
        @param articles: Iterable of L{Article} objects.
        """
        with open(self.artlist_path, "ab") as artlist:
            # Add new articles to article list and database
            for art in articles:
                pmid = str(art.pmid)
                if pmid not in self.artdb: 
                    artlist.write("%s %d\n" % (
                        pmid, Date2Integer(art.date_completed)))
                    self.artdb[pmid] = art
            self.artdb.sync()
