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
from path import path
import sys

from mscanner.medline.Databases import ArticleData
from mscanner.medline.Updater import Updater
from mscanner.configuration import rc
from mscanner.core import iofuncs

def regenerate():
    """Recreate feature map, feature stream, and featuredatabase"""
    updater = Updater.Defaults()
    updater.regenerate()


def add_directory():
    """Add new articles to MScanner databases by parsing XML files in
    a Medline directory."""
    logging.info("Updating MScanner from " + rc.medline.relpath())
    updater = Updater.Defaults()
    updater.add_directory(rc.medline, save_delay=3)


def add_pickle(pickle):    
    """Add articles to MScanner database from a gzip'd pickle.
    
    @param pickle: Path to a gzip'd pickled list of Article objects."""
    pickle = path(pickle)
    logging.debug("Updating database from pickle %s", pickle.basename())
    with closing(GzipFile(pickle, "rb")) as zf:
        articles = cPickle.load(zf)
    updater = Updater.Defaults()
    updater.add_articles(articles)

    
def make_pickle(filename):
    """Write a Gzip Pickle of Article objects for later use with L{add_pickle}.
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
            art = adata.artdb[str(pmid)]
            articles.append(art)
        except KeyError:
            logging.debug("%d not found in database", pmid)
    logging.debug("Writing articles to pickle %s", pickle.basename())
    with closing(GzipFile(pickle, "wb")) as zf:
        cPickle.dump(articles, zf, protocol=2)


if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger(logfile=False)
    locals()[sys.argv[1]](*sys.argv[2:])
