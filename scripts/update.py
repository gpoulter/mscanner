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

from mscanner.medline.Updater import Updater
from mscanner.configuration import rc
from mscanner.core import iofuncs

featurespaces = [
    "feats_mesh_qual_issn",
    #"feats_wmqia",
    #"feats_wmqia_filt",
    #"feats_iedb_word", 
    #"feats_iedb_concat",
    #"feats_word_fold",
    #"feats_word_num", 
    #"feats_word_strip",
]


def regenerate():
    """Regenerate the FeatureMap, FeatureStream, FeatureVectors
    and article list - but only those that have been deleted from dist."""
    updater = Updater.Defaults(featurespaces)
    for fdata in updater.fdata_list:
        fdata.regenerate(updater.artdb)
    return updater


def medline():
    """Add new articles to MScanner databases by parsing XML files in
    a Medline directory."""
    logging.info("Updating MScanner from " + rc.medline.relpath())
    updater = Updater.Defaults(featurespaces)
    updater.add_directory(rc.medline, save_delay=0)
    return updater


def vacuum(mincount):
    """Remove low-frequency features from the index, remapping the feature vectors.
    @param mincount: Keep features with at least this many occurrences.
    """
    mincount = int(mincount)
    logging.info("Vacuuming features with counts less than %d", mincount)
    updater = Updater.Defaults(featurespaces)
    for fdata in updater.fdata_list:
        fdata.vacuum(mincount=mincount)
    return updater


def load_pickles(*pickles):    
    """Add articles to MScanner database from gzip'd pickles, each
    of which contains a list of Article objects."""
    updater = Updater.Defaults(featurespaces)
    for pickle in pickles:
        pickle = path(pickle)
        logging.debug("Updating database from pickle %s", pickle.basename())
        with closing(GzipFile(pickle, "rb")) as zf:
            articles = cPickle.load(zf)
        updater.add_articles(articles)
    return updater

    
def save_pickle(filename):
    """Write a Gzip Pickle of Article objects for later use with L{load_pickles}.
    If filename is x.txt, the output will be in x.pickle.gz.
    
    @param filename: Path to list of PubMed IDs.
    """
    from mscanner.medline import Shelf
    filename = path(filename)
    pickle = filename.stripext() + ".pickle.gz"
    logging.debug("Retrieving articles for file %s", filename.basename())
    artdb = Shelf.open(rc.articles_home/rc.articledb,'r')
    articles = []
    for pmid in iofuncs.read_pmids(filename):
        try:
            articles.append(artdb[str(pmid)])
        except KeyError:
            logging.debug("%d not found in database", pmid)
    logging.debug("Writing articles to pickle %s", pickle.basename())
    with closing(GzipFile(pickle, "wb")) as zf:
        cPickle.dump(articles, zf, protocol=2)


if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger()
    locals()[sys.argv[1]](*sys.argv[2:])
