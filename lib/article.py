"""Handle Articles, processed files, features, feature occurrence counts

@author: Graham Poulter
                                     

StatusFile -- One way posting of program status to a file
Article -- Stores attributes of an article
FileTracker -- Track processed files (to avoid re-parsing), with on-disk persistence
FeatureMapping -- Mapping between string features and 16-bit feature IDs

countFeatures() -- Get feature occurrence counts over a set of documents
getArticles() -- Retrieve Articles from a DB, caching results in a Pickle
readPMIDs() -- Get PMIDs (and scores too) from a text file
makeBackup() -- Create a .recovery file
removeBackup() -- Remove a .recovery file


"""

from array import array
from codecs import open
import cPickle
import dbshelve
from itertools import izip

import logging as log
import os
import time

from path import path
import numpy

def countFeatures(nfeats, featdb, docids):
    """Givent number of features,

    @param nfeats: Number of distinct features (size of array)

    @param featdb: Mapping from document ID to list of feature IDs

    @param docids: List of document IDs whose features are to be counted

    @return: Array of length nfeats, with occurrence count of each feature
    """
    import numpy
    counts = numpy.zeros(nfeats, dtype=numpy.int32)
    for docid in docids:
        counts[featdb[docid]] += 1
    return counts

def readPMIDs(filename, include=None, exclude=None, withscores=False):
    """Read PubMed IDs one per line from filename.

    @param filename: Path to file containing one PubMed ID per line,
    with optional score after the PubMed ID. File format ignores blank
    lines and lines starting with #, and only parses the line up to
    the first whitespace character.

    @param include: Only return members of this set

    @param exclude: Do not return members of this set
    
    @param withscores: Also read the score after the PubMed ID on each line.
    
    @returns: An iterator over PubMed ID, or (PubMed ID, Score) if withscores==True
    """
    if not isinstance(filename,path) or not filename.exists():
        raise ValueError("File %s does not exist" % filename)
    count = 0
    for line in file(filename, "r"):
        sline = line.strip()
        if sline == "" or sline.startswith("#"):
            continue
        splits = sline.split()
        pmid = int(splits[0])
        if include is not None and pmid not in include:
            fname, ext = filename.splitext()
            (fname+".broken"+ext).write_lines([str(pmid)],append=True)
            log.warn("PMID %d is not a member of include" % pmid)
            continue
        if exclude is not None and pmid in exclude:
            fname, ext = filename.splitext()
            (fname+".exclude"+ext).write_lines([str(pmid)],append=True)
            log.warn("PMID %d is a member of exclude" % pmid)
            continue
        if withscores:
            yield pmid, float(splits[1])
        else:
            yield pmid
        count += 1
    if count == 0:
        raise RuntimeError("Did not succeed in reading any PMIDs from %s" % filename)

def writePMIDScores(filename, pairs):
    """Write PubMed IDs and their scores to a file
    
    @pairs: Iteratable over (PubMed ID, Score) pairs
    """
    oscores = sorted(pairs, key=lambda x:x[1], reverse=True)
    filename.write_lines("%d %f" % x for x in oscores)

def runMailer(smtp_server, mailer):
    """Read e-mail addresses from mail file and send them a
    message saying MScanner is available.

    @param mailer: Path object to file of e-mail addresses
    """
    if not mailer.isfile():
        return
    import smtplib
    import logging
    server = smtplib.SMTP(smtp_server)
    server.set_debuglevel(0)
    fromaddr = "nobody@mscanner.stanford.edu"
    for email in mailer.lines(retain=False):
        logging.debug("Sending availability alert to %s", email)
        msg = "From: %s\r\nTo: %s\r\n\r\n" % (fromaddr, email)
        msg += """
Someone - hopefully you - requested that an alert be sent when
MScanner (http://mscanner.stanford.edu) completed its current request.

MScanner has completed its request and is currently available for
 query or validation operations at the time of sending of this e-mail.

(MScanner classifies Medline citations by training a Naive Bayes
classifier on their metadata, making it useful for information
retrieval using example citations)
"""
        try:
            server.sendmail(fromaddr, email, msg)
        except Exception:
            logging.debug("Failed to send to %s", email)
    server.quit()
    mailer.remove()

def getArticles(article_db_path, pmidlist_path):
    """Get Article objects given a file of PubMed IDs.

    @param article_db_path: Path to a berkeley DB mapping PubMed IDs
    to Article objects.

    @param pmidlist_path: Path to a text file listing one PubMed ID
    per line.

    @return: List of Article objects in the order given in the text
    file.

    @note: The first time it is called for a given pmidlist_path, the
    results are cached in a .pickle appended to pmidlist_path, and
    later calls simply use the cached results.
    """
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        return cPickle.load(file(cache_path, "rb"))
    pmids = readPMIDFile(pmidlist_path)
    artdb = dbshelve.open(article_db_path, 'r')
    articles = [artdb[str(p)] for p in pmids]
    cPickle.dump(articles, file(cache_path,"wb"), protocol=2)
    artdb.close()
    return articles

def makeBackup(filepath):
    """Make a .recovery backup of a file.

    @raise RuntimeError: if the backup already exists.
    """
    if filepath == "None": return
    oldfile = filepath+".recover"
    if oldfile.exists():
        raise RuntimeError("Recovery required: %s exists" % oldfile)
    if filepath.exists():
        filepath.rename(oldfile)

def removeBackup(filepath):
    """Remove the .recovery backup of a file.

    @raise RuntimeError: if the backup does not exist.
    """
    if filepath == "None": return
    if not filepath.exists():
        raise RuntimeError("Recovery required: %s does not exist (backup removal requested)" % filepath)
    oldfile = filepath+".recover"
    if oldfile.exists():
        oldfile.remove()

class StatusFile:
    """Class to manage a status file, for one way posting of program
    status to listeners."""

    def __init__(self, filename, dataset, total=0):
        """Create status file with PID, progress, total, dataset and start time"""
        if filename is None:
            raise RuntimeError("No status file given")
        if filename.exists():
            raise RuntimeError("Status file %s already exists" % str(filename))
        self.filename = filename
        self.dataset = dataset
        self.progress = 0
        self.total = total
        self.start = time.time()
        self.write()

    def __del__(self):
        """Remove status file on deletion"""
        if hasattr(self,"filename") and self.filename.isfile():
            self.filename.remove()

    def __str__(self):
        """Return status file contents"""
        return "%d\n%d\n%d\n%s\n%s\n" % (
            os.getpid(), self.progress, self.total, self.dataset, str(self.start))

    def write(self):
        """Write status file to disk"""
        self.filename.write_text(str(self))
        
    def update(self, progress, total=None):
        """Update progress of status file.  If progress is None, set to total."""
        if total is not None:
            self.total = total
        if progress is None:
            self.progress = self.total
        else:
            self.progress = progress
        log.info("Completed %d out of %d", self.progress, self.total)
        self.write()

    __call__ = update
    
class FileTracker(set):
    """Class which tracks processed files.

    It accepts all paths, but membership is checked according to basename()
    """

    def __init__(self, trackfile=None):
        """Initialise, specifying path to write the list of files"""
        self.trackfile = trackfile
        if isinstance(trackfile,path) and trackfile.isfile():
            self.update(trackfile.lines(retain=False))

    def dump(self):
        """Write the list of tracked files, one per line"""
        if self.trackfile is None:
            return
        makeBackup(self.trackfile)
        self.trackfile.write_lines(sorted(self))
        removeBackup(self.trackfile)

    def add(self, fname):
        """Add a file to the tracker"""
        set.add(self, fname.basename())

    def toprocess(self, files):
        """Given a list of files return sorted list of those not in the tracker"""
        return sorted(f for f in files if f.basename() not in self)

class FeatureMapping:
    """Curates a database of string features, providing methods to map
    between a feature string and an integer feature ID.

    @note: The type of a feature could be 'mesh', 'qual', 'issn',
    'year', or 'author'.  A given feature string could have more than
    one type.

    @ivar feats: List mapping ID to (feature string, feature type)
    @ivar counts: List mapping ID to number of occurrences
    @ivar feat2id: {type:{feature:id}} mapping
    """

    def __init__(self, featfile=None):
        """Initialise the database

        @param featfile: Path to text file with list of terms
        """
        self.featfile = featfile
        self.feats = []
        self.feat2id = {}
        self.counts = []
        if featfile is not None and self.featfile.exists():
            self.load()
        
    def load(self):
        """Load feature mapping mapping from file

        @note: file format is '%(feature)\t%(type)\t%(count)\n'
        """
        self.feats = []
        self.feat2id = {}
        f = open(self.featfile, "rb", "utf-8")
        for fid, line in enumerate(f):
            feat, ftype, count = line.strip().split("\t")
            self.feats.append((feat,ftype))
            self.counts.append(int(count))
            if ftype not in self.feat2id:
                self.feat2id[ftype] = {feat:fid}
            else:
                self.feat2id[ftype][feat] = fid
        f.close()

    def dump(self):
        """Write the feature mapping to disk"""
        makeBackup(self.featfile)
        f = open(self.featfile, "wb", "utf-8")
        for (feat, ftype), count in zip(self.feats, self.counts):
            f.write(feat+"\t"+ftype+"\t"+str(count)+"\n")
        f.close()
        removeBackup(self.featfile)

    def __getitem__(self, featid):
        """Return (feature string, feature type) given feature ID"""
        return self.feats[featid]

    def __len__(self):
        """Return number of distinct features"""
        return len(self.feats)

    def getFeatures(self, feature_ids):
        """Return (string,type) for features given feature id list"""
        return [self.feats[fid] for fid in feature_ids]

    def featureCounts(self):
        """Return array with number of occurrences of each feature"""
        return numpy.array(self.counts, numpy.int32)

    def featureTypeMask(self, exclude_types):
        """Get a mask for excluded features

        @param exclude_types: Types of features to exclude

        @return: Boolean array for excluded features (but returns None if exclude_types is None)
        """
        if exclude_types is None:
            return None
        exclude_feats = numpy.zeros(len(self.feats), dtype=numpy.bool)
        for ftype in exclude_types:
            for fid in self.feat2id[ftype].itervalues():
                exclude_feats[fid] = True
        return exclude_feats

    def getFeatureIds(self, features, ftype, count=False):
        """Get feature IDs corresponding to feature strings.

        @note: Dynamically creates new features IDs and types as necessary.

        @param features: List of feature strings to convert

        @param ftype: Type of the feature strings

        @param count: Whether or not to add to occurrence counts
        
        @return: List of feature IDs
        """
        result = []
        if ftype not in self.feat2id:
            self.feat2id[ftype] = {}
        fdict = self.feat2id[ftype]
        if count:
            for feat in features:
                if feat not in fdict:
                    featid = len(self.feats)
                    self.feats.append((feat,ftype))
                    self.counts.append(1)
                    fdict[feat] = featid
                else:
                    self.counts[fdict[feat]] += 1
                result.append(fdict[feat])
        else:
            for feat in features:
                if feat not in fdict:
                    raise ValueError("Unknown feature in read-only feature mapping: " + feat)
                else:
                    result.append(fdict[feat])
        return result

class Article:
    """A simple wrapper for parsed Medline articles.
    
    Instance variables are same as the constructor parameters.

    """
    def __init__(self,
                 pmid=0,
                 title="",
                 abstract="",
                 journal="",
                 issn="",
                 year=0,
                 meshterms=None,
                 authors=None):
        """
        Initialise a new article with the given parameters.
        
        @param pmid: Integer PubMed ID or MEDLINE UI of the article.
        @param title: Title of the article
        @param abstract: Abstract of the article
        @param journal: Medline abbreviated journal title
        @param issn: ISSN code for the journal
        @param year: Year of publication
        @param meshterms: Set of Mesh terms associated with article.
        @param authors: Set of (initials,lastname) pairs of article authors
        """
        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract
        self.journal = journal
        self.issn = issn
        self.year = year 
        self.meshterms = meshterms
        if meshterms is None:
            self.meshterms = list()
        self.authors = authors
        if authors is None:
            self.authors = list()

    def __repr__(self):
        import pprint
        pp = pprint.PrettyPrinter()
        astr = "Article(pmid=%d,\ntitle=%s,\nabstract=%s,\njournal=%s\nissn=%s\nyear=%s\nmeshterms=%s\nauthors=%s)\n"
        return astr % (
            self.pmid,
            repr(self.title),
            repr(self.abstract),
            repr(self.journal),
            repr(self.issn),
            repr(self.year),
            pp.pformat(self.meshterms),
            pp.pformat(self.authors))

