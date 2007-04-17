"""For updating the databases of articles and features (feature
mapping, db of article->features and stream of article->features).

                                   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html

http://www.gnu.org/copyleft/gpl.html
"""

from bsddb import db
import cPickle
import gzip
import logging as log
import numpy as nx
import os
from path import path
import xml.etree.cElementTree as ET

from article import Article
import dbshelve
from featuredb import FeatureDatabase, FeatureStream
from featuremap import FeatureMapping
from utils import FileTracker

def parse(source):
    """Convert XML file into Article objects
    
    @param source: File-like object containing XML and MedlineCitation elements
    
    @return: Iterator over Article objects
    """
    context = ET.iterparse(source, events=("start", "end"))
    context = iter(context)
    event, root = context.next()
    for event, record in context:
        if event == "end" and record.tag == "MedlineCitation":
            if record.get("Status") == "MEDLINE":
                result = Article()
                result.pmid = int(record.findtext("PMID"))
                art = record.find("Article")
                result.issn = art.findtext("Journal/ISSN")
                result.year = art.findtext("Journal/JournalIssue/PubDate/Year")
                if result.year is not None:
                    result.year = int(result.year)
                result.title = art.findtext("ArticleTitle")
                result.abstract = art.findtext("Abstract/AbstractText")
                result.authors = [(a.findtext("Initials"), a.findtext("LastName")) 
                                  for a in art.findall("AuthorList/Author")]
                result.journal = record.findtext("MedlineJournalInfo/MedlineTA")
                for heading in record.findall("MeshHeadingList/MeshHeading"):
                    descriptor = heading.findtext("DescriptorName")
                    quals = [ q.text for q in heading.findall("QualifierName") ]
                    result.meshterms.append(tuple([descriptor] + quals))
                yield result
            root.clear()

class MedlineCache:
    """Manages the database of medline abstracts.  Function is to
    update the cache with new articles, and retrieve articles from the
    cache.
    """

    def __init__(
        self, featmap, db_env_home,
        article_db, feature_db, feature_stream, article_list,
        processed_path, narticles_path, use_transactions=True):
        """Initialse a cache of the results of parsing medline.

        @param featmap: A FeatureMapping object for mapping string features to IDs
        @param db_env_home: Path to DB home directory 
        @param article_db: Path to article database
        @param feature_db: Path to feature database
        @param feature_stream: Path to feature stream file
        @param article_list: Path to list of article PMIDs
        @param narticles_path: Path to file containing the total number of PMIDs
        @param processed_path: Path to list of processed files
        @param use_transactions: If false, disable transaction engine
        """
        self.db_env_home = db_env_home
        self.featmap = featmap
        self.article_db = article_db
        self.feature_db = feature_db
        self.feature_stream = feature_stream
        self.article_list = article_list
        self.processed_path = processed_path
        self.narticles_path = narticles_path
        self.use_transactions = use_transactions
        self.recover = False

    def makeDBEnv(self):
        """Initialise DB environment for transactions"""
        if not self.db_env_home.isdir():
            self.db_env_home.mkdir()
        dbenv = db.DBEnv()
        dbenv.set_lg_max(128*1024*1024) # 128Mb log files
        dbenv.set_tx_max(1) # 1 transaction at a time
        dbenv.set_cachesize(0, 8*1024*1024) # 8Mb shared cache
        flags = db.DB_INIT_MPOOL|db.DB_CREATE
        if self.use_transactions:
            flags |= db.DB_INIT_TXN
        if self.recover:
            flags |= db.DB_RECOVER # might use db.DB_RECOVER_FATAL
        dbenv.open(self.db_env_home, flags)
        return dbenv

    def putArticles(self, articles, dbenv):
        """Store Article objects and feature lists the databases

        @param articles: Iterable of Article objects
        @param dbenv: Database environment to use
        """
        log.info("Starting transaction to add articles")
        txn = None
        if self.use_transactions:
            txn = dbenv.txn_begin()
        try:
            artdb = dbshelve.open(self.article_db, dbenv=dbenv, txn=txn)
            meshfeatdb = FeatureDatabase(self.feature_db, dbenv=dbenv, txn=txn)
            featstream = FeatureStream(file(self.feature_stream,"ab"))
            if not self.narticles_path.isfile():
                narticles = len(meshfeatdb)
            else:
                narticles = int(self.narticles_path.text())
            pmidlist = []
            for art in articles:
                # Refuse to add duplicates
                if art.pmid in meshfeatdb: 
                    continue
                # Store article, adding it to list of documents
                artdb[str(art.pmid)] = art
                pmidlist.append(str(art.pmid))
                # Get MeSH headings, qualifiers and ISSN from article
                headings = list()
                quals = list()
                issns = list()
                for term in art.meshterms:
                    headings.append(term[0])
                    if(len(term)>1):
                        for q in term[1:]:
                            if q not in quals:
                                quals.append(q)
                if art.issn:
                    issns = [art.issn]
                # Add features to feature mapping
                featids = self.featmap.addArticle(mesh=headings, qual=quals, issn=issns)
                meshfeatdb.setitem(art.pmid, featids, txn)
                featstream.write(art.pmid, featids)
            artdb.close()
            meshfeatdb.close()
            featstream.close()
            self.article_list.write_lines(pmidlist, append=True)
            self.narticles_path.write_text(str(narticles+len(pmidlist))+"\n")
            self.featmap.dump()
            if txn is not None:
                txn.commit()
        except Exception, e:
            if txn is not None:
                log.error("Aborting Transaction: Error %s", e)
                txn.abort()
            raise
        else:
            if txn is not None:
                log.info("Committed transaction")
            return len(pmidlist)
            
    def updateCacheFromDir(self, medlinedir, save_delay=5):
        """Updates the cache given that medlinedir contains .xml.gz
        file to add to the cache and that we should save the inverse
        document after processing each savesteps files."""
        import time
        filenames = medlinedir.files("*.xml") + medlinedir.files("*.xml.gz")
        tracker = FileTracker(self.processed_path)
        toprocess = tracker.toprocess(filenames)
        dbenv = self.makeDBEnv()
        for idx, filename in enumerate(toprocess):
            log.info("Adding to cache: file %d out of %d (%s)", idx+1, len(toprocess), filename.name)
            for t in xrange(save_delay):
                log.debug("Saving in %d seconds...", save_delay-t)
                time.sleep(1)
            log.debug("Parsing XML file %s", filename.basename())
            if filename.endswith(".gz"):
                infile = gzip.open(filename, 'r')
            else:
                infile = file(filename, 'r')
            numadded = self.putArticles(parse(infile), dbenv)
            log.debug("Added %d articles", numadded)
            tracker.add(filename)
            tracker.dump()
            log.info("Completed file %d out of %d (%s)", idx+1, len(toprocess), filename.name)
        dbenv.close()
