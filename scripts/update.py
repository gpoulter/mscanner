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

import cPickle
import logging
import sys

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
    updater.add_directory(rc.medline, save_delay=0)


def add_pickle(pickle):    
    """Add articles to MScanner database from a pickle.
    
    @param pickle: Path to a pickled list of Article objects."""
    logging.info("Updating MScanner from " + pickle)
    with open(pickle , "rb") as f:
        articles = cPickle.load(f)
    updater = Updater.Defaults()
    updater.add_articles(articles)


if __name__ == "__main__":
    # Call the named function with provided arguments
    iofuncs.start_logger(logfile=False)
    locals()[sys.argv[1]](*sys.argv[2:])