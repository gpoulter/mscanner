#!/usr/bin/env python

"""Update the MScanner database with new articles

Usage::
    python update.py [somepickle]

If a path to a pickle is given, load Article objects from the Pickle into
the MScanner database.

With no arguments, look for new XML files in the Medline path and add their
contents to the database. 
"""

                                     
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
import logging as log
import sys

from mscanner import medline, featuremap
from mscanner.configuration import rc, initLogger


def update_mscanner(pickle=None):
    """Look for MScanner updates in the rc.medline direcotry
    
    @param pickle: Path to a pickle of Article objects to be added instead.
    """
    medcache = medline.MedlineCache(
            featuremap.FeatureMapping(rc.featuremap),
            rc.db_env_home,
            rc.articledb,
            rc.featuredb,
            rc.featurestream,
            rc.articlelist,
            rc.processed,
            rc.narticles,
            rc.use_transactions,
            )
    # Load articles from a pickle
    if pickle is not None:
        log.info("Updating MScanner from " + pickle )
        articles = cPickle.load(open(pickle , "rb"))
        dbenv = medcache.create_dbenv()
        medcache.add_articles(articles, dbenv)
        dbenv.close()
    # Parse articles from XML directory
    else:
        log.info("Updating MScanner from " + rc.medline.relpath())
        medcache.add_directory(rc.medline, rc.save_delay)


if __name__ == "__main__":
    initLogger(logfile=False)
    if len(sys.argv) == 1:
        update_mscanner()
    elif len(sys.argv) == 2:
        update_mscanner(sys.argv[1])
    else:
        print __doc__
