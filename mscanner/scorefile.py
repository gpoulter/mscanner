"""Reads/writes files of PubMed IDs/scores, but also loads articles
from database and prints an error page for when no PubMed IDs are valid."""

from __future__ import with_statement
from __future__ import division

                                     
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

import logging as log
from contextlib import closing
import numpy as nx

from mscanner.configuration import rc
from mscanner import featuredb, featuremap
from mscanner.support import dbshelve


class Databases:
    """Connections to the MScanner databases

    The environment needs to be reloaded when databases are updated,
    because L{featmap} and L{article_list} will have changed on disk.

    @ivar featdb: Mapping from PubMed ID to list of features

    @ivar featmap: L{FeatureMapping} between feature names and IDs
    
    @ivar artdb: Mmapping from PubMed ID to Article object
    """
    
    def __init__(self):
        """Constructor for setting attributes to be used by the remaining
        methods."""
        log.info("Loading article databases")
        self.featdb = featuredb.FeatureDatabase(rc.featuredb, 'r')
        self.featmap = featuremap.FeatureMapping(rc.featuremap)
        self.artdb = dbshelve.open(rc.articledb, 'r')


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



def read_pmids(filename, 
              include=None, 
              exclude=None, 
              broken_name=None, 
              exclude_name=None, 
              withscores=False):
    """Read PubMed IDs one per line from filename.

    @param filename: Path to file with PubMed IDs, formatted one PubMed ID per
    line, with optional score after the PubMed ID. Blank lines and lines
    starting with # are ignored, as is data after the PMID and score.
    
    @param include: Only return members of this set (other PubMed IDs are
    considered "broken").

    @param broken_name: File to write non-included ("broken") PubMed IDs
    
    @param exclude: Do not return members of this set
    
    @param exclude_name: File to write excluded PubMed IDs
    
    @param withscores: Also read the score after the PubMed ID on each line.
    
    @returns: Iterator over PubMed ID, or (Score, PubMed ID) if withscores
    True. """
    count = 0
    broken = []
    excluded = []
    if filename.isfile(): 
        with open(filename, "r") as infile:
            for line in infile:
                sline = line.strip()
                if sline == "" or sline.startswith("#"):
                    continue
                splits = sline.split()
                pmid = int(splits[0])
                if include is not None and pmid not in include:
                    broken.append(pmid)
                    continue
                if exclude is not None and pmid in exclude:
                    excluded.append(pmid)
                    continue
                if withscores:
                    yield float(splits[1]), pmid
                else:
                    yield pmid
                count += 1
    if broken_name is not None:
        with open(broken_name, "w") as f:
            f.write("\n".join(str(s) for s in broken))
    if exclude_name is not None:
        with open(exclude_name, "w") as f:
            f.write("\n".join(str(s) for s in excluded))
    log.debug("Got %d PubMed IDs from %s", count, filename.basename())

    
def write_scores(filename, pairs):
    """Write scores and PubMed IDs to file, in decreasing order of score.
    @param pairs: Iterable of (score, PMID)     
    """
    from path import path
    path(filename).write_lines(
        "%-10d %f" % (p,s) for s,p in sorted(pairs, reverse=True))
    
    
def load_articles(article_db_path, pmidlist_path):
    """Retrieve Article objects given a file of PubMed IDs.

    @param article_db_path: Path to a berkeley DB mapping PubMed IDs
    to Article objects.

    @param pmidlist_path: Path to a text file listing one PubMed ID per line.

    @return: List of Article objects in the order given in the text file.

    @note: The first called with a given PMID list caches the results in a
    .pickle, and later calls load the pickle."""
    import cPickle
    from mscanner.support import dbshelve
    from path import path # used in the line below
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        with open(cache_path, "rb") as f:
            return cPickle.load(f)
    pmids = read_pmids(pmidlist_path)
    with closing(dbshelve.open(article_db_path, "r")) as artdb:
        articles = [artdb[str(p)] for p in pmids]
    with open(cache_path, "wb") as f:
        cPickle.dump(articles, f, protocol=2)
    return articles


def no_valid_pmids_page(filename, pmids):
    """Print an error page when no valid PMIDs were found
    @param filename: Path to output file
    @param pmids: List of any provided PMIDs (all invalid)
    """
    import logging
    from Cheetah.Template import Template
    from mscanner.support import utils
    from mscanner.configuration import rc
    logging.warning("No valid PubMed IDs were found!")
    with utils.FileTransaction(filename, "w") as ft:
        page = Template(file=str(rc.templates/"notfound.tmpl"))
        page.notfound_pmids = pmids
        page.respond(ft)
