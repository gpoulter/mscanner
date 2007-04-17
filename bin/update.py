#!/usr/bin/env python

"""Update MEDLINE database with articles

Usage: update.py [somepickle]

Without arguments: Parse .xml[.gz] files in configured Medline path
and add Article objects to databse.

With one argument: Load Article objects from a Pickle and them to
database.

                                   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import cPickle
import logging as log
import sys

from mscanner import article
from mscanner import configuration as c
from mscanner.featuremap import FeatureMapping
from mscanner.utils import StatusFile, runMailer
from mscanner.medline import MedlineCache

if __name__ == "__main__":
    # Usage information
    if len(sys.argv)>1:
        if "-h" in sys.argv or "--help" in sys.argv:
            print __doc__
            sys.exit(0)
    
    # Initialise database
    medcache = MedlineCache(
            FeatureMapping(c.featuremap),
            c.db_home,
            c.articledb,
            c.featuredb,
            c.featurestream,
            c.articlelist,
            c.processed,
            c.narticles,
            c.use_transactions,
            )
    
    statfile = StatusFile(c.statfile, dataset="UPDATING DATABASE", start=c.timestamp)
    try: 
        # Load articles from a pickle
        if len(sys.argv) == 2:
            print "Loading articles from " + sys.argv[1]
            articles = cPickle.load(file(sys.argv[1], "rb"))
            dbenv = medcache.makeDBEnv()
            medcache.putArticleList(articles, dbenv)
            dbenv.close()
            
        # Parse articles from XML directory
        else:
            log.info("Starting update from %s" % c.medline.relpath())
            medcache.updateCacheFromDir(c.medline, c.save_delay)
    finally:
        log.debug("Cleaning up")
        del statfile
        runMailer(c.smtp_server, c.emails_path)
