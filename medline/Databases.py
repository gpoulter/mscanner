"""For consumers of the database, this opens L{FeatureDatabase}, 
L{FeatureMapping} and the article list"""

from __future__ import with_statement
import logging as log
import numpy as nx
from contextlib import closing
from mscanner.configuration import rc
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.medline.FeatureDatabase import FeatureDatabase
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


class Databases:
    """Connections to the MScanner databases

    The environment needs to be reloaded when databases are updated,
    because L{featmap} and L{article_list} will have changed on disk.

    @ivar featdb: Mapping from PubMed ID to list of features

    @ivar featmap: L{FeatureMapping} between feature names and IDs
    (by default, featmap[id] == feature string)
    
    @ivar artdb: Mmapping from PubMed ID to Article object
    """
    
    def __init__(self):
        """Constructor for setting attributes to be used by the remaining
        methods."""
        log.info("Loading article databases")
        self.featdb = FeatureDatabase(rc.featuredb, 'r')
        self.featmap = FeatureMapping(rc.featuremap)
        self.artdb = Shelf.open(rc.articledb, 'r')


    @property
    def article_list(self):
        """Array with the PubMed IDs in the database 
        
        Array may be over 16 million members long, so the property and will
        take a while to load the first time it is accessed."""
        try:
            return self._article_list
        except AttributeError: 
            log.info("Loading article list")
            self._article_list = nx.array([int(x) for x in rc.articlelist.lines()])
            return self._article_list


    def close(self):
        """Closes the feature and article databases"""
        self.featdb.close()
        self.artdb.close()
    __del__ = close


def load_articles(db_path, pmids_path):
    """Retrieve Articles listed in a file of PubMed IDs

    @param db_path: Path to a berkeley DB mapping PubMed IDs
    to Article objects.

    @param pmids_path: Path to a text file listing one PubMed ID per line.

    @return: List of Article objects in the order given in the text file.

    @note: The first called with a given PMID list caches the results in a
    .pickle, and later calls load the pickle."""
    import cPickle
    from mscanner.medline import Shelf
    from contextlib import closing
    from mscanner.core.iofuncs import read_pmids
    from path import path # used in the line below
    cache_path = path(pmids_path + ".pickle")
    if cache_path.isfile():
        with open(cache_path, "rb") as f:
            return cPickle.load(f)
    pmids = read_pmids(pmids_path)
    with closing(Shelf.open(db_path, "r")) as artdb:
        articles = [artdb[str(p)] for p in pmids]
    with open(cache_path, "wb") as f:
        cPickle.dump(articles, f, protocol=2)
    return articles