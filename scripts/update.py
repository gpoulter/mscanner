#!/usr/bin/env python

"""Update the MScanner database with new articles

Usage::
    python update.py [somepickle]

If a path to a pickle is given, load Article objects from the Pickle into
the MScanner database.

With no arguments, look for new XML files in the Medline path and add their
contents to the database. 
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


def update_dir():
    """Add articles to MScanner databases by parsing XML files in
    a Medline directory."""
    updater = Updater.Defaults()
    logging.info("Updating MScanner from " + rc.medline.relpath())
    updater.add_directory(rc.medline, save_delay=0)


def update_pickle(pickle):    
    """Add articles to MScanner database from a pickle.
    @param pickle: Path to a pickled list of Article objects."""
    logging.info("Updating MScanner from " + pickle )
    with open(pickle , "rb") as f:
        articles = cPickle.load(f)
    updater = Updater.Defaults()
    updater.add_articles(articles)


if __name__ == "__main__":
    iofuncs.start_logger(logfile=False)
    if len(sys.argv) == 1:
        update_dir()
    elif len(sys.argv) == 2:
        update_pickle(sys.argv[1])
    else:
        print __doc__
    logging.shutdown()
