#!env python

"""Update MEDLINE database with articles

@author: Graham Poulter
                                   

Usage: update.py [somepickle]

No arguments: Parse .xml[.gz] files in configuration.medline_path and
add Article's to databse.

One argument: Load Article objects from a Pickle and add to database

"""

# Builtin modules
import cPickle
import logging as log
import sys
# Custom modules
import configuration as c
import article
import medline
import xmlparse

# Usage information
if len(sys.argv)>1:
    if "-h" in sys.argv or "--help" in sys.argv:
        print __doc__
        sys.exit(0)

# Initialise database
medcache = medline.MedlineCache(
        article.FeatureMapping(c.featuremap),
        xmlparse.ArticleParser(),
        c.db_home,
        c.articledb,
        c.featuredb,
        c.articlelist,
        c.processed,
        c.use_transactions,
        )

statfile = article.StatusFile(c.statfile, "UPDATING DATABASE")
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
    article.runMailer(c.smtp_server, c.mailer)
    
        
