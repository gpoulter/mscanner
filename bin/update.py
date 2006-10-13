#!/usr/bin/env python

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
from configuration import update as u, medline as m
from xmlparse import ArticleParser
from medline import MedlineCache
from article import FeatureMapping

# Usage information
if len(sys.argv)>1:
    if "-h" in sys.argv or "--help" in sys.argv:
        print __doc__
        sys.exit(0)

# Initialise database
parser = ArticleParser(u.mesh_synonyms, u.mesh_excludes)
meshdb = FeatureMapping(m.meshdb)
med = MedlineCache(
        meshdb,
        parser,
        m.db_home,
        m.articledb,
        m.featuredb,
        m.termcounts,
        m.processed,
        )
    
# Load Articles from a Pickle
if len(sys.argv) == 2:
    print "Loading articles from " + sys.argv[1]
    articles = cPickle.load( file( sys.argv[1], "rb" ) )
    exclude = cPickle.load( file( u.mesh_excludes, "rb" ) )
    for article in articles:
        article.meshterms -= exclude
    dbenv = med.makeDBEnv()
    med.putArticleList(articles, dbenv)
    dbenv.close()

# Parse articles from XML directory
else:
    log.info("Starting update from %s" % u.medline.relpath())
    med.updateCacheFromDir( u.medline, u.save_delay )
