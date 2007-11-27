"""Syncrhonising around the article and feature databases. They allow the
database to be loaded at the same time, and for articles are added
simultaneously to all of them databases."""

from __future__ import with_statement
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
    """
    
    def __init__(self, featmap, featdb, fstream, ftype, dbenv=None):
        """Constructor. 
        @param dbenv: Optional Berkeley database environment"""
        self.featmap = FeatureMapping(featmap, ftype=ftype)
        self.featuredb = FeatureDatabase(featdb, dbenv=dbenv, ftype=ftype)
        self.fstream = FeatureStream(fstream, ftype=ftype)


    @staticmethod
    def Defaults_MeSH(dbenv=None):
        """Construct using the RC defaults for MeSH feature space."""
        return FeatureData(rc.featuremap_mesh, rc.featuredb_mesh, 
                           rc.featurestream_mesh, nx.uint16, dbenv)


    @staticmethod
    def Defaults_All(dbenv=None):
        """Construct using RC defaults for the MeSH+Abstract feature space."""
        return FeatureData(rc.featuremap_all, rc.featuredb_all, 
                           rc.featurestream_all, nx.uint32, dbenv)


    def close(self):
        """Shut down the databases"""
        self.featmap.dump()
        self.featuredb.close()
        self.fstream.close()


    def add_articles(self, articles):
        """Add articles to the databases and feature maps. 
        @param articles: Iterable of (PMID, date, features), as (string,
        YYYYMMDD integer, dictionary). The features dict maps string feature
        type to list of feature strings."""
        for pmid, date, features in articles:
            if pmid not in self.featuredb:
                featvec = self.featmap.add_article(features)
                self.featuredb[pmid] = featvec
                self.fstream.additem(pmid, date, featvec)
        self.featmap.dump()
        self.featuredb.sync()
        self.fstream.flush()



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
            self._article_list = nx.array(
                [int(x.split()[0]) for x in self.artlist_path.lines()])
            return self._article_list


    def add_articles(self, articles):
        """Add articles to the database and article list
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


    def load_articles(self, pmids_path):
        """Return Article objects given a file listing PubMed IDs, caching
        the results in a pickle.  The pickle has a ".pickle" on top of
        L{pmids_path}.
        
        @param pmids_path: Path to a text file with one PubMed ID per line.
        
        @return: List of Article objects (same order as text file).
        """
        import cPickle
        from path import path
        from mscanner.core.iofuncs import read_pmids
        cache_path = path(pmids_path + ".pickle")
        if cache_path.isfile():
            with open(cache_path, "rb") as f:
                return cPickle.load(f)
        else:
            articles = [self.artdb[str(p)] for p in read_pmids(pmids_path)]
            with open(cache_path, "wb") as f:
                cPickle.dump(articles, f, protocol=2)
            return articles