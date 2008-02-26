#!/usr/bin/env python

"""Update or regenerate the MScanner databases.

Usage::
    ./update.py function_name arg1 arg2 [...]

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

from contextlib import closing
import cPickle
from gzip import GzipFile
import logging
import numpy as nx
from path import path
import sys

from mscanner.medline.ArticleData import ArticleData
from mscanner.medline.Updater import Updater
from mscanner.configuration import rc
from mscanner.core import iofuncs

featurespaces = [
    ("feats_mesh_qual_issn", nx.uint16), 
    ("feats_wmqia", nx.uint32),
    ("feats_wmqia_filt", nx.uint32),
    ("feats_iedb_word", nx.uint32),
    ("feats_iedb_concat", nx.uint32),
    ("feats_word_nofold", nx.uint32),
    ("feats_word_num", nx.uint32),
    ("feats_word_strip", nx.uint32),
]


def regenerate():
    """Regenerate the FeatureMap, FeatureStream, FeatureDatabase
    and article list - but only those that have been deleted from dist."""
    updater = Updater.Defaults(featurespaces)
    updater.regenerate()


def medline():
    """Add new articles to MScanner databases by parsing XML files in
    a Medline directory."""
    logging.info("Updating MScanner from " + rc.medline.relpath())
    updater = Updater.Defaults(featurespaces)
    updater.add_directory(rc.medline, save_delay=3)


def load_pickles(*pickles):    
    """Add articles to MScanner database from gzip'd pickles, each
    of which contains a list of Article objects."""
    updater = Updater.Defaults(featurespaces)
    for pickle in pickles:
        pickle = path(pickle)
        logging.debug("Updating database from pickle %s", pickle.basename())
        with closing(GzipFile(pickle, "rb")) as zf:
            articles = cPickle.load(zf)
        updater.add_articles(articles, sync=False)
    updater.sync() # Only syncing writes FeatureMapping to disk.
    updater.close() # Closing won't write the FeatureMapping.

    
def save_pickle(filename):
    """Write a Gzip Pickle of Article objects for later use with L{load_pickles}.
    If filename is x.txt, the output will be in x.pickle.gz.
    
    @param filename: Path to list of PubMed IDs.
    """
    filename = path(filename)
    pickle = filename.stripext() + ".pickle.gz"
    logging.debug("Retrieving articles for file %s", filename.basename())
    adata = ArticleData.Defaults()
    articles = []
    for pmid in iofuncs.read_pmids(filename):
        try:
            articles.append(adata.artdb[str(pmid)])
        except KeyError:
            logging.debug("%d not found in database", pmid)
    logging.debug("Writing articles to pickle %s", pickle.basename())
    with closing(GzipFile(pickle, "wb")) as zf:
        cPickle.dump(articles, zf, protocol=2)


if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger(logfile=False)
    locals()[sys.argv[1]](*sys.argv[2:])
