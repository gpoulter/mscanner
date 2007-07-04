"""For updating the databases of articles and features (feature
mapping, db of article->features and stream of article->features).

parse() -- Turn Medline XML
MedlineCache -- For consistent updating of the database/features
FileTracker -- Track processed files (to avoid re-parsing), with persistence

                                   
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

from bsddb import db
import gzip
import logging as log
from path import path
import xml.etree.cElementTree as ET

from mscanner.support import dbshelve
from mscanner.article import Article
from mscanner.featuredb import FeatureDatabase, FeatureStream
from mscanner.featuremap import FeatureMapping

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
    """Class for updating the Article DB, FeatureMapping,
    FeatureDatabase, FeatureStream, PMID list, and FileTracker.
    
    @note: Databases are opened and closed for each putArticles() transaction
    to limit the impact of interruptions (how do I make it atomic?)
    
    Passed through constructor:
        @ivar featmap: A FeatureMapping object for mapping string features to IDs
        @ivar db_env_home: Path to DB home directory 
        @ivar article_db: Path to article database
        @ivar feature_db: Path to feature database
        @ivar feature_stream: Path to feature stream file
        @ivar article_list: Path to list of article PMIDs
        @ivar narticles_path: Path to file containing the total number of PMIDs
        @ivar processed_path: Path to list of processed files
        @ivar use_transactions: If false, disable transaction engine
    """

    def __init__(
        self, 
        featmap, 
        db_env_home,
        article_db, 
        feature_db, 
        feature_stream, 
        article_list,
        processed_path, 
        narticles_path, 
        use_transactions=True):
        """Initialse the database
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
        
class FileTracker(set):
    """A persistent set for tracking of processed files.
    
    @note: Overloads "add" to only add basenames, and facilities to 
    load from and dump to file.

    It accepts all paths, but membership is checked according to basename()
    """

    def __init__(self, trackfile=None):
        """Initialise, specifying path to write the list of files"""
        self.trackfile = trackfile
        self.trackfile_new = trackfile+".new"
        if isinstance(trackfile,path) and trackfile.isfile():
            self.update(trackfile.lines(retain=False))

    def dump(self):
        """Write the list of tracked files, one per line"""
        if self.trackfile is None:
            return
        self.trackfile_new.write_lines(sorted(self))
        self.trackfile_new.rename(self.trackfile)

    def add(self, fname):
        """Add a file to the tracker"""
        set.add(self, fname.basename())

    def toprocess(self, files):
        """Given a list of files return sorted list of those not in the tracker"""
        return sorted(f for f in files if f.basename() not in self)


