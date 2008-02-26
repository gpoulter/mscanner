"""Wraps files with the Article database (maps PubMed ID to Article), article
list (of PubMed IDs) and file with number of articles in the databse"""

from __future__ import with_statement
from contextlib import closing
import logging
import numpy as nx

from mscanner.configuration import rc
from mscanner.medline import Shelf 
from mscanner.medline.FeatureStream import DateAsInteger
from mscanner import endpath

                                     
__author__ = "Graham Poulter"                                        
__license__ = """GPThis program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


class ArticleData:
    """Wraps the files containing Article data - namely the Article database,
    and the list of PubMed IDs and dates.  
    
    @ivar artdb: Mapping from PubMed ID to Article object.
    
    @ivar artlist_path: Path to file listing PubMed ID and record date.
    """
    
    def __init__(self, artdb_path, artlist_path, artcount_path):
        """Constructor for setting attributes to be used by the remaining
        methods."""
        self.artdb = Shelf.open(artdb_path)
        self.artlist_path = artlist_path
        self.artcount_path = artcount_path


    @staticmethod
    def Defaults():
        """Instantiate L{ArticleData} using default paths"""
        base = rc.articles_home
        # Create new article repository if necessary
        if not base.exists(): base.makedirs()
        return ArticleData(base/rc.articledb, 
                           base/rc.articlelist,
                           base/rc.articlecount)


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
            logging.info("Loading article list from %s", endpath(self.artlist_path))
            if self.artlist_path.exists():
                # Save memory by not constructing a list() object
                self._article_list = nx.zeros(self.article_count, nx.uint32)
                with open(self.artlist_path, "rb") as f:
                    for idx, line in enumerate(f):
                        self._article_list[idx] = int(line.split()[0])
            else:
                self._article_list = nx.array([], nx.uint32)
            return self._article_list
        
        
    @property
    def article_count(self):
        """Number of PubMed IDs in the database. We read these from
        L{artcount_path} because len(self.artdb) is very slow."""
        try:
            return self._article_count
        except AttributeError:
            logging.info("Loading article count from  %s", endpath(self.artcount_path))
            if self.artcount_path.exists():
                with open(self.artcount_path, "rb") as f:
                    self._article_count = int(f.readline())
            else:
                self._article_count = 0
            return self._article_count


    def regenerate_artlist(self):
        """Regenerate the list of articles (and number of articles) from the
        article database. Flushes to disk at the end."""
        if self.artlist_path.exists():
            logging.info("Not regenerating article list/count: the file already exists.")
            return
        if hasattr(self, "_article_list"):
            del self._article_list
        if hasattr(self, "_article_count"):
            del self._article_count
        with open(self.artlist_path, "wb") as f:
            self._article_count = 0
            for pmid, art in self.artdb.iteritems():
                self._article_count += 1 
                if self._article_count % 10000 == 0:
                    logging.debug("Wrote %d to the article list." % numdocs)
                f.write("%s %08d\n" % (pmid, DateAsInteger(art.date_completed)))
            f.flush()
        with open(self.artcount_path, "wb") as f:
            f.write("%d\n" % self._article_count)


    def add_articles(self, articles):
        """Add new articles to both the database and the article list.  Flushes
        to disk at the end.
        @param articles: Iterable of L{Article} objects.
        """
        # Load article count
        self.article_count
        # New articles get added to list, count and database
        with open(self.artlist_path, "ab") as artlist:
            for art in articles:
                pmid = str(art.pmid)
                if pmid not in self.artdb: 
                    artlist.write("%s %d\n" % (
                        pmid, DateAsInteger(art.date_completed)))
                    self.artdb[pmid] = art
                    self._article_count += 1
            self.artdb.sync()
            artlist.flush()
        # Write new article count
        with open(self.artcount_path, "wb") as f:
            f.write("%d\n" % self._article_count)
