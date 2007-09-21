#!/usr/bin/env python

"""Update the MScanner database with new articles

Usage::
    python update.py [somepickle]

Without no arguments, parses .xml[.gz] files in configured Medline path and add
Article objects to the database.

With one argument it loads Article objects from a Pickle and adds them to the
database."""

                                     
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

from mscanner import article, medline, featuremap
from mscanner.configuration import rc, initLogger

def main():
    # Initialise database
    initLogger(logfile=False)
    log.info("Initialising databases")
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
    if len(sys.argv) == 2:
        log.info("Loading articles from " + sys.argv[1])
        articles = cPickle.load(open(sys.argv[1], "rb"))
        dbenv = medcache.create_dbenv()
        medcache.add_articles(articles, dbenv)
        dbenv.close()
    # Parse articles from XML directory
    else:
        log.info("Starting update from %s" % rc.medline.relpath())
        medcache.add_directory(rc.medline, rc.save_delay)
    log.debug("Cleaning up")

if __name__ == "__main__":
    # Usage information
    if len(sys.argv)>1:
        if "-h" in sys.argv or "--help" in sys.argv:
            print __doc__
            sys.exit(0)
    main()
