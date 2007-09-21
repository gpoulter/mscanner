"""Performs reading/writing of various file types

Functions in this file:
    - L{read_pmids}: Get PMIDs (and scores too) from a text file
    - L{write_scores}: Write PMIDs and scores to a text file
    - L{load_articles}: Retrieve Articles from a DB, caching results in a Pickle
    - L{read_descriptor}: Read paramaters from a descriptor file
    - L{write_descriptor}: Write paramaters to a descriptor file
"""

from __future__ import with_statement
                                     
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
        from mscanner.configuration import rc
        from mscanner import featuredb, featuremap
        from mscanner.support import dbshelve
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
    from path import path
    from mscanner.support import dbshelve
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        return cPickle.load(open(cache_path, "rb"))
    pmids = read_pmids(pmidlist_path)
    artdb = dbshelve.open(article_db_path, 'r')
    articles = [artdb[str(p)] for p in pmids]
    cPickle.dump(articles, open(cache_path,"wb"), protocol=2)
    artdb.close()
    return articles


def parsebool(s):
    """Handler for converting strings to booleans"""
    if isinstance(s, basestring):
        s = s.strip()
        if s == "0" or s == "False":
            return False
        elif s == "1" or s == "True":
            return True
        else:
            raise ValueError("Failed to parse boolean: %s" % s)
    else:
        return bool(s)


descriptor_keys = dict(
    alpha=float,
    code=str,
    dataset=str,
    delcode=str,
    hidden=parsebool,
    limit=int,
    nfolds=int,
    numnegs=int,
    operation=str,
    threshold=float,
    starttime=float,
    timestamp=float,)


def read_descriptor(fpath):
    """Reads a descriptor file, returning a dictionary of parameters.

    Each descriptor line is formatted as "#key = value". Reading stops at the
    first line that does not start with '#'. Valid keys are in
    L{descriptor_keys}. The same file can be used with L{read_pmids}, which
    ignores lines beginning with '#'.

    @return: Storage object, with additional "_filename" key that contains
    fpath."""
    f = open(fpath, "r")
    line = f.readline()
    from mscanner.support.storage import Storage
    result = Storage()
    while line.startswith("#"):
        key, value = line[1:].split(" = ",1)
        value = descriptor_keys[key](value.strip())
        result[key] = value
        line = f.readline()
    result["_filename"] = fpath
    f.close()
    return result


def write_descriptor(fpath, pmids, params):
    """Write parameters and PubMed IDs to the descriptor file.
    
    @param fpath: File to write

    @param pmids: List of PubMed IDs, may be None
    
    @param params: Key-value dictionary to write. Values are converted with
    str(). Only keys that have a member in descriptor_keys are actually written
    to the file."""
    f = open(fpath, "w")
    for key, value in params.iteritems():
        if key in descriptor_keys: 
            f.write("#" + key + " = " + str(value) + "\n")
    if pmids is not None:
        for pmid in pmids:
            f.write(str(pmid)+"\n")
    f.close()


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
        page.rc = rc
        page.notfound_pmids = pmids
        page.respond(ft)
