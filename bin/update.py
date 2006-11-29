#!env python

"""Update MEDLINE database with articles

@author: Graham Poulter
                                   

Usage: update.py [somepickle]

No arguments: Parse .xml[.gz] files in configuration.medline_path and
add Article's to databse.

One argument: Load Article objects from a Pickle and add to database

"""

import cPickle
import logging as log
import sys

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
parser = xmlparse.ArticleParser(c.mesh_synonyms, c.mesh_excludes)
featmap = article.FeatureMapping(c.featuremap)
medcache = medline.MedlineCache(
        featmap,
        parser,
        c.db_home,
        c.articledb,
        c.featuredb,
        c.articlelist,
        c.termcounts,
        c.processed,
        c.use_transactions,
        )
    
# Load Articles from a Pickle
if len(sys.argv) == 2:
    print "Loading articles from " + sys.argv[1]
    articles = cPickle.load(file(sys.argv[1], "rb"))
    exclude = cPickle.load(file(c.mesh_excludes, "rb"))
    for article in articles:
        article.meshterms -= exclude
    dbenv = medcache.makeDBEnv()
    medcache.putArticleList(articles, dbenv)
    dbenv.close()

# Parse articles from XML directory
else:
    log.info("Starting update from %s" % c.medline.relpath())
    medcache.updateCacheFromDir(c.medline, c.save_delay)
