"""Wraps the FeatureMapping, FeatureVectors and FeatureStream objects,
so that they are opened and closed together, and articles are added
to all three representations at the same time."""

from __future__ import with_statement
from contextlib import closing
import logging
import numpy as nx

from mscanner.configuration import rc
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.medline.FeatureDatabase import FeatureDatabase
from mscanner.medline.FeatureStream import FeatureStream, DateAsInteger
from mscanner import endpath


                                     
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
    """Wraps the L{FeatureMapping}, L{FeatureDatabase} and L{FeatureStream}
    objects. Different L{FeatureData} instances are possible for the same
    L{ArticleData} collection, depending on the choice of feature extraction
    method.

    @ivar featmap: L{FeatureMapping} between feature names and feature IDs.
    
    @ivar featuredb: Mapping from PubMed ID to list of features.
    
    @ivar fstream: L{FeatureStream} of (PMID, date, feature vector).
    
    @ivar featurespace: Name of a method of L{Article} which will generate
    the features, which is also the subdirectory in rc.articles_home where
    the feature databases are to be kep.
    
    @ivar rdonly: If True, open all databases read-only.
    """
    
    def __init__(self, numtype, featmap, featdb, fstream, featurespace, rdonly=True):
        """Constructor.
        @param numtype: Numpy type for FeatureDatabase (uint16 or uint32).
        """
        self.rdonly = rdonly
        logging.debug("Loading feature mapping from %s", endpath(featmap))
        self.featmap = FeatureMapping(featmap)
        self.featuredb = FeatureDatabase(numtype, featdb, "r" if rdonly else "c")
        self.fstream = FeatureStream(fstream, rdonly)
        self.featurespace = featurespace


    @staticmethod 
    def Defaults(featurespace, numtype, rdonly=True):
        """Initialise L{FeatureData} using standard file paths. numtype,
        L{featurespace} and L{rdonly} are as in the constructor."""
        base = rc.articles_home / featurespace
        # Create feature database directory if necessary
        if not base.exists(): base.makedirs()
        return FeatureData(numtype, base/rc.featuremap, base/rc.featuredb, 
                           base/rc.featurestream, featurespace, rdonly)
    

    def close(self):
        """Shut down the databases"""
        self.featuredb.close()
        self.fstream.close()
        self.featmap.close()


    def add_articles(self, articles):
        """Incrementally add new articles to the existing feature 
        database, stream and feature map.  
        
        @param articles: Iterable Article objects."""
        if self.rdonly:
            raise NotImplementedError("Attempt to write to read-only databases")
        logging.debug("Adding articles to %s", endpath(self.featmap.filename.dirname()))
        for article in articles:
            pmid = article.pmid
            if pmid not in self.featuredb:
                date = DateAsInteger(article.date_completed)
                features = getattr(article, self.featurespace)()
                self.featmap.add_article(features)
                featvec = self.featmap.make_vector(features)
                self.featuredb[str(pmid)] = featvec
                self.fstream.additem(pmid, date, featvec)
        self.featuredb.sync()
        self.fstream.flush()
        self.featmap.con.commit()
        
        
    def regenerate(self, articles):
        """Regenerate feature map, feature stream and feature database, but
        only if they have been deleted. FeatureMapping is flushed to disk only
        at the end.
        
        @param articles: Iterable of all Article objects in the L{ArticleData}."""
        if self.rdonly:
            raise NotImplementedError("Failed: may not write read-only databases.")
        do_featmap = len(self.featmap) == 1
        do_stream = self.fstream.filename.size == 0
        do_featuredb = len(self.featuredb) == 0
        if not (do_featmap or do_stream or do_featuredb):
            logging.info("Done: nothing to do as databases already have data.")
            return
        if do_featmap:
            logging.info("Regen feature map %s.", endpath(self.featmap.filename))
            if not (do_stream and do_featuredb):
                raise ValueError("Cannot regenerate feature map without regenerating feature stream and database")
        if do_stream: 
            logging.info("Regen feature stream %s.", endpath(self.fstream.filename))
        if do_featuredb: 
            logging.info("Regen feature database %s.", endpath(self.featuredb.filename))
        for count, article in enumerate(articles):
            pmid = article.pmid
            date = DateAsInteger(article.date_completed)
            features = getattr(article, self.featurespace)()
            if count % 10000 == 0:
                logging.debug("Regenerated %d articles.", count)
                self.featmap.con.commit()
            if do_featmap:
                self.featmap.add_article(features)
            if do_stream or do_featuredb:
                featvec = self.featmap.make_vector(features)
                if do_stream:
                    self.fstream.additem(pmid, date, featvec)                
                if do_featuredb:
                    self.featuredb[str(pmid)] = featvec
        self.featuredb.sync()
        self.fstream.flush()
        self.featmap.con.commit()
