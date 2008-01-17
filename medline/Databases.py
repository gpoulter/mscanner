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
from mscanner.medline.FeatureStream import FeatureStream, DateAsInteger
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


def _endof(path):
    """Return last two parts of the path, for displaying paths when the full
    path is too long, but just the filename looks ambiguous."""
    return "/".join(path.splitall()[-2:])


class FeatureData:
    """Wraps the files handling feature ID information - namely the
    L{FeatureMapping}, L{FeatureDatabase} and L{FeatureStream}.

    @ivar featmap: L{FeatureMapping} between feature names and feature IDs.
    
    @ivar featuredb: Mapping from PubMed ID to list of features.
    
    @ivar fstream: L{FeatureStream} of (PMID, date, feature vector).
    
    @ivar extractor: Function taking an L{Article} parameter and returning
    feature dictionary as {space:[features]}.
    
    @ivar featurespace: Name of a method of L{Article} which will generate
    the features, which is also the subdirectory in rc.articles_home where
    the feature databases are to be kep.
    
    @ivar ftype: Numpy type for FeatureDatabase (uint16 or uint32).
        
    @ivar rdonly: If True, open all databases read-only.
    """
    
    def __init__(self, ftype, featmap, featdb, fstream, featurespace, rdonly=True):
        self.ftype = ftype
        self.rdonly = rdonly
        logging.debug("Loading feature mapping from %s", _endof(featmap))
        self.featmap = FeatureMapping(featmap)
        self.featuredb = FeatureDatabase(ftype, featdb, "r" if rdonly else "c")
        self.fstream = FeatureStream(fstream, rdonly)
        self.featurespace = featurespace


    @staticmethod 
    def Defaults(featurespace, ftype, rdonly=True):
        """Initialise L{FeatureData} using standard file paths. L{ftype},
        L{featurespace} and L{rdonly} are as in the constructor."""
        base = rc.articles_home / featurespace
        # Create feature database directory if necessary
        if not base.exists(): base.makedirs()
        return FeatureData(ftype, base/rc.featuremap, base/rc.featuredb, 
                           base/rc.featurestream, featurespace, rdonly)
    

    def close(self):
        """Shut down the databases"""
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
        
        @param articles: Iterable Article objects."""
        if self.rdonly:
            raise NotImplementedError("Attempt to write to read-only databases")
        for article in articles:
            pmid = article.pmid
            if pmid not in self.featuredb:
                date = DateAsInteger(article.date_completed)
                features = getattr(article, self.featurespace)()
                self.featmap.add_article(features)
                featvec = self.featmap.get_vector(features)
                self.featuredb[str(pmid)] = featvec
                self.fstream.additem(pmid, date, featvec)
        
        
    def regenerate(self, articles):
        """Regenerate feature map, feature stream and feature database, but
        only if they have been deleted. FeatureMapping is flushed to disk only
        at the end.
        
        @param articles: Iterable of all Article objects in the L{ArticleData}."""
        if self.rdonly:
            raise NotImplementedError("Failed: may not write read-only databases.")
        do_featmap = len(self.featmap) == 0
        do_stream = self.fstream.filename.size == 0
        do_featuredb = len(self.featuredb) == 0
        if not (do_featmap or do_stream or do_featuredb):
            logging.info("Done: nothing to do as databases already have data.")
            return
        if do_featmap: 
            logging.info("Regen feature map %s.", _endof(self.featmap.filename))
        if do_stream: 
            logging.info("Regen feature stream %s.", _endof(self.fstream.filename))
        if do_featuredb: 
            logging.info("Regen feature database %s.", _endof(self.featuredb.filename))
        for count, article in enumerate(articles):
            pmid = article.pmid
            date = DateAsInteger(article.date_completed)
            features = getattr(article, self.featurespace)()
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
    
    def __init__(self, artdb_path, artlist_path):
        """Constructor for setting attributes to be used by the remaining
        methods."""
        self.artdb = Shelf.open(artdb_path)
        self.artlist_path = artlist_path


    @staticmethod
    def Defaults():
        """Instantiate L{ArticleData} using default paths"""
        base = rc.articles_home
        # Create new article repository if necessary
        if not base.exists(): base.makedirs()
        return ArticleData(base/rc.articledb, base/rc.articlelist)


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
            logging.info("Loading article list from %s", _endof(self.artlist_path))
            if self.artlist_path.exists():
                with open(self.artlist_path, "rb") as f:
                    self._article_list = nx.array(
                        [int(t.split()[0]) for t in f], nx.uint32)
            else:
                self._article_list = nx.array([], nx.uint32)
            return self._article_list


    def regenerate_artlist(self):
        """Regenerate the list of articles from the article database. Flushes
        to disk at the end."""
        if self.artlist_path.exists():
            logging.info("Not regenerating article list: the file already exists.")
            return
        if hasattr(self, "_article_list"):
            del self._article_list
        with open(self.artlist_path, "wb") as f:
            for idx, (pmid, art) in enumerate(self.artdb.iteritems()):
                if idx % 10000 == 0:
                    logging.debug("Wrote %d to the article list." % idx)
                f.write("%s %08d\n" % (pmid, DateAsInteger(art.date_completed)))
            f.flush()


    def add_articles(self, articles):
        """Add new articles to both the database and the article list.  Flushes
        to disk at the end.
        @param articles: Iterable of L{Article} objects.
        """
        with open(self.artlist_path, "ab") as artlist:
            # Add new articles to article list and database
            for art in articles:
                pmid = str(art.pmid)
                if pmid not in self.artdb: 
                    artlist.write("%s %d\n" % (
                        pmid, DateAsInteger(art.date_completed)))
                    self.artdb[pmid] = art
            self.artdb.sync()
            artlist.flush()
