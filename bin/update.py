#!/usr/bin/env python

"""Update MEDLINE database with articles

Usage: update.py [somepickle]

Without arguments: Parse .xml[.gz] files in configured Medline path
and add Article objects to databse.

With one argument: Load Article objects from a Pickle and them to
database.

                                   
"""

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

from mscanner import article
from mscanner.configuration import rc, initLogger
from mscanner.featuremap import FeatureMapping
from mscanner.medline import MedlineCache

if __name__ == "__main__":
    # Usage information
    if len(sys.argv)>1:
        if "-h" in sys.argv or "--help" in sys.argv:
            print __doc__
            sys.exit(0)
    # Initialise database
    initLogger(logfile=False)
    log.info("Initialising databases")
    medcache = MedlineCache(
            FeatureMapping(rc.featuremap),
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
        articles = cPickle.load(file(sys.argv[1], "rb"))
        dbenv = medcache.makeDBEnv()
        medcache.putArticleList(articles, dbenv)
        dbenv.close()
    # Parse articles from XML directory
    else:
        log.info("Starting update from %s" % rc.medline.relpath())
        medcache.updateCacheFromDir(rc.medline, rc.save_delay)
    log.debug("Cleaning up")
