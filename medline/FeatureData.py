"""Wraps the FeatureMapping, FeatureVectors and FeatureStream objects,
so that they are opened and closed together, and articles are added
to all three representations at the same time."""

from __future__ import with_statement
from contextlib import closing
import logging
import numpy as nx
import sys

from mscanner.configuration import rc
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.medline.FeatureVectors import FeatureVectors
from mscanner.medline.FeatureStream import FeatureStream, DateAsInteger, vb_encode
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
    """Wraps the L{FeatureMapping}, L{FeatureVectors} and L{FeatureStream}
    objects. There may be multiple L{FeatureData} indexes for the same database
    of articles, depending on the choice of feature extraction method.

    @ivar featmap: L{FeatureMapping} between feature names and feature IDs.
    
    @ivar featuredb: L{FeatureVectors} to look up feature vectors by PubMed ID.
    
    @ivar fstream: L{FeatureStream} of (PMID, date, feature vector) records.
    
    @ivar featurespace: Name of a method of L{Article} which will generate
    the features, which is also the subdirectory in rc.articles_home where
    the feature databases are to be kep.
    
    @ivar rdonly: If True, treat all databases as read-only.
    """
    
    def __init__(self, featmap, featdb, fstream, featurespace, rdonly=True):
        """Constructor.
        @param featmap: Path to FeatureMapping.
        @param featdb: Path to FeatureVectors.
        @param fstream: Path to FeatureStream.
        """
        logging.debug("Loading features from %s", endpath(featmap.dirname()))
        self.rdonly = rdonly
        self.featmap = FeatureMapping(featmap)
        self.featuredb = FeatureVectors(featdb)
        self.fstream = FeatureStream(fstream, rdonly)
        self.featurespace = featurespace


    @staticmethod 
    def Defaults(featurespace, rdonly=True):
        """Initialise L{FeatureData} using standard file paths. L{featurespace}
        and L{rdonly} are as in the constructor."""
        # Create index directory if necessary
        base = rc.articles_home / featurespace
        if not base.exists(): base.makedirs()
        return FeatureData(base/rc.featuremap, base/rc.featuredb, 
                           base/rc.featurestream, featurespace, rdonly)
    

    def close(self):
        """Shut down the databases"""
        self.featuredb.close()
        self.fstream.close()
        self.featmap.close()


    def commit(self):
        """Commit transactions"""
        self.featmap.con.commit()
        self.featuredb.con.commit()
        self.fstream.flush()


    def add_articles(self, articles):
        """Incrementally add new articles to the existing feature 
        database, stream and feature map.  
        
        @param articles: Iterator over Article objects."""
        if self.rdonly:
            raise NotImplementedError("Attempt to write to read-only databases")
        logging.debug("Adding articles to %s", endpath(self.featmap.filename.dirname()))
        for article in counter(articles):
            pmid = article.pmid
            if pmid not in self.featuredb:
                date = DateAsInteger(article.date_completed)
                features = getattr(article, self.featurespace)()
                featvec = vb_encode(self.featmap.add_article(features))
                self.featuredb.add_record(pmid, date, featvec)
                self.fstream.additem(pmid, date, featvec)
        self.commit()


    def regenerate(self, artdb):
        """Regenerate feature map, feature stream and feature database, but
        only if they have been deleted. FeatureMapping is flushed to disk only
        at the end.
        
        @param artdb: Dictionary of Article objects keyed by PubMed ID."""
        if self.rdonly:
            raise NotImplementedError("Failed: may not write read-only databases.")
        do_featmap = len(self.featmap) == 1
        do_stream = self.fstream.filename.size == 0
        do_featuredb = len(self.featuredb) == 0
        if not (do_featmap or do_stream or do_featuredb):
            logging.info("Regen: nothing to do as databases already have data.")
            return
        # Regenerate map, db, stream
        if do_featmap:
            logging.info("Regenerating map,db,stream %s.", endpath(self.featmap.filename.dirname()))
            if not (do_stream and do_featuredb):
                raise ValueError("Cannot regenerate feature map without doing stream/database as well.")
            self.add_articles(artdb.itervalues())
        # Regenerate feature stream from database
        elif do_stream: 
            logging.info("Regenerating FeatureStream %s.", endpath(self.fstream.filename))
            for pmid, date, featvec in counter(self.featuredb.iteritems(decode=False)):
                self.fstream.additem(pmid, date, featvec)
        # Regenerate feature database from feature stream
        elif do_featuredb: 
            logging.info("Regenerating FeatureVectors %s.", endpath(self.featuredb.filename))
            for pmid, date, featvec in counter(self.fstream.iteritems(decode=False)):
                self.featuredb.add_record(pmid, date, featvec)
        self.commit()
        logging.info("Index regeneration complete.")


def counter(iter):
    """Passes through the iterator, printing '.' on stdout after every 1000
    items.  Adds a newline at the end."""
    count = 0
    for x in iter:
        yield x
        count += 1
        if count >= 1000:
            sys.stdout.write(".")
            sys.stdout.flush()
            count = 0
    sys.stdout.write("\n")
    sys.stdout.flush()
