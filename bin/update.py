#!env python

"""Update MEDLINE database with articles

@author: Graham Poulter
                                   

Usage: update.py [somepickle]

Without arguments: Parse .xml[.gz] files in configured Medline path and
add Article's to databse.

With one argument: Load Article objects from a Pickle and them to database.

"""

import cPickle
import logging as log
import sys

import configuration as c
from featuremap import FeatureMapping
from utils import StatusFile, runMailer
from medline import MedlineCache

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
        c.use_transactions,
        )

statfile = StatusFile(c.statfile, "UPDATING DATABASE")
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
    del statfile
    runMailer(c.smtp_server, c.mailer)
    
        
