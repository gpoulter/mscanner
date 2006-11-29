"""Handle Articles, processed files, features, feature occurrence counts

@author: Graham Poulter
                                     

chooseRandomLines() -- Get random lines from a file
readPMIDFile() -- Get PMIDs reliably from a text file
getArticles() -- Retrieve Articles from a DB, caching results in a Pickle
makeBackup() -- Create a .recovery file
removeBackup() -- Remove a .recovery file

StatusFile -- One way posting of program status to a file
Article -- Stores attributes of an article
FileTracker -- Track processed files (to avoid re-parsing), with on-disk persistence
FeatureMapping -- Mapping between string features and 16-bit feature IDs
TermCounts -- Store number of occurrences of each feature, + total occurrences and documents

"""

from array import array
import cPickle
import dbshelve
from itertools import izip
import logging as log
import os
from path import path
import time

def getArticles(article_db_path, pmidlist_path):
    """Return list of Article's, given the path to an article database
    and the path to a file with a list of pubmed IDs.

    The first time it is called for a given pmidlist_path, the results
    are cached in a .pickle appended to pmidlist_path, and later calls
    simply use the cached results.
    """
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        return cPickle.load(file(cache_path, "rb"))
    pmids = readPMIDFile(pmidlist_path)
    artdb = dbshelve.open(article_db_path, 'r')
    articles = [ artdb[str(p)] for p in pmids ]
    cPickle.dump(articles, file(cache_path,"wb"), protocol=2)
    artdb.close()
    return articles

def makeBackup(filepath):
    """Make a .recovery backup of a file.

    Raise RuntimeError if the backup already exists.
    """
    if filepath == "None": return
    oldfile = filepath+".recover"
    if oldfile.exists():
        raise RuntimeError("Recovery required: %s exists" % oldfile)
    if filepath.exists():
        filepath.rename(oldfile)

def removeBackup(filepath):
    """Remove the .recovery backup of a file.

    Raise RuntimeError if the backup does not exist.
    """
    if filepath == "None": return
    if not filepath.exists():
        raise RuntimeError("Recovery required: %s does not exist (backup removal requested)" % filepath)
    oldfile = filepath+".recover"
    if oldfile.exists():
        oldfile.remove()

def chooseRandomLines(infile_path, outfile_path, N):
    """Choose N random lines from infile_path and print them to outfile_path"""
    from random import randint
    lines = file(infile_path,"r").readlines()
    size = len(lines)
    outfile = file(outfile_path,"w")
    if N > size:
        raise ValueError("N > length of file")
    i = 0
    while i < N:
        i += 1
        r = randint(0,size-1)
        outfile.write(lines[r])
        lines[r] = lines[size-1]
        size = size - 1

def readPMIDFile(filename, allpmids=None):
    """Read PubMed IDs one per line from filename.

    File format ignores blank lines and lines starting with #, and only
    parses the line up to the first whitespace character.

    If allpmids is provided, we only return those PMID integers that
    satisfy 'pmid in allpmids'.
    """
    if not isinstance(filename,path) or not filename.exists():
        raise ValueError("File %s does not exist" % filename)
    count = 0
    for line in file(filename,"r"):
        if line.strip() != "" and not line.startswith("#"):
            pmid = int(line.split()[0])
            if allpmids is None or pmid in allpmids:
                yield pmid
                count += 1
            else:
                log.warn("Failed to find PMID %d in allpmids" % pmid)
    if count == 0:
        raise RuntimeError("Did not succeed in reading any PMIDs from %s" % filename)

class StatusFile:
    """Class to manage a status file, for one way posting of program
    status to listeners."""

    def __init__(self, filename, dataset, total=0):
        """Create status file with PID, progress, total, dataset and start time"""
        self.filename = filename
        self.dataset = dataset
        self.progress = 0
        self.total = total
        self.start = time.time()
        if self.filename.exists():
            raise RuntimeError("Status file already exists")
        self.write()

    def __del__(self):
        """Remove status file on deletion"""
        self.filename.remove()

    def __str__(self):
        """Return status file contents"""
        return "%d\n%d\n%d\n%s\n%s\n" % (
            os.getpid(),
            self.progress,
            self.total,
            self.dataset,
            str(self.start))

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
        self.write()

class Article:
    """A simple wrapper for parsed Medline articles

    @ivar pmid: Integer PubMed ID or MEDLINE UI of the article.
    @ivar title: Title of the article
    @ivar abstract: Abstract of the article
    @ivar journal: ISO abbreviated journal name
    @ivar issn: ISSN code for the journal
    @ivar year: Year of publication
    @ivar meshterms: Set of Mesh terms associated with article.
    @ivar authors: Set of (initials,lastname) pairs of article authors
    @ivar chemicals: Set of chemicals associated with article
    """
    def __init__(self,
                 pmid=0,
                 title="",
                 abstract="",
                 journal="",
                 issn="",
                 year=0,
                 meshterms=None,
                 authors=None,
                 chemicals=None):
        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract
        self.journal = journal
        self.issn = issn
        self.year = year 
        self.meshterms = meshterms
        if meshterms is None:
            self.meshterms = set()
        self.authors = authors
        if authors is None:
            self.authors = list()
        self.chemicals = chemicals
        if chemicals is None:
            self.chemicals = set()

    def __repr__(self):
        import pprint
        pp = pprint.PrettyPrinter()
        astr = "Article(pmid=%d,\ntitle=%s,\nabstract=%s,\njournal=%s\nissn=%s\nyear=%s\nmeshterms=%s\nauthors=%s\nchemicals=%s)\n"
        return astr % (self.pmid,
                       repr(self.title),
                       repr(self.abstract),
                       repr(self.journal),
                       repr(self.issn),
                       repr(self.year),
                       pp.pformat(self.meshterms),
                       pp.pformat(self.authors),
                       pp.pformat(self.chemicals),)

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

    @note: The type of a feature could be 'mesh', 'issn', 'year', or
     'author'.  A given feature string could have more than one type.

    @ivar feats: List mapping ID to (feature string, feature type)
    @ivar feat2id: {type:{feature:id}} mapping
    """

    def __init__(self, featfile=None):
        """Initialise the database

        @param featfile: Path to text file with list of terms
        """
        self.featfile = featfile
        self.feats = []
        self.feat2id = {}
        self.load()

    def load(self):
        """Load the term mapping from disk"""
        self.feats = []
        self.types = []
        self.feat2id = {}
        if not isinstance(self.featfile,path) or not self.featfile.exists():
            return
        f = file(self.featfile, "r")
        for fid, pair in enumerate(f):
            pair = pair.strip().decode("utf-8").split("\t")
            if len(pair) == 1:
                pair.append(None)
            feat, ftype = pair
            self.feats.append((feat,ftype))
            if ftype not in self.feat2id:
                self.feat2id[ftype] = {feat:fid}
            else:
                self.feat2id[ftype][feat] = fid
        f.close()

    def dump(self):
        """Write the feature mapping to disk"""
        if self.featfile is None:
            return
        makeBackup(self.featfile)
        f = self.featfile.open("wb",)
        for feat, ftype in self.feats:
            f.write(feat.encode("utf-8"))
            if ftype is not None:
                f.write("\t"+ftype)
            f.write("\n")
        f.close()
        removeBackup(self.featfile)

    def __getitem__(self, featid):
        """Return feature string given feature ID"""
        return self.feats[featid]

    def __len__(self):
        """Return number of distinct features"""
        return len(self.feats)

    def getFeatures(self, featids):
        """Return (string,type) for features given feature id list"""
        return [ self.feats[fid] for fid in featids ]

    def getFeatureIds(self, feats, ftype=None):
        """Get term IDs corresponding to the given list of
        terms. Dynamically creates new features IDs as necessary.

        @param feats: List of features strings to convert

        @param ftype: Type of the features being added
        
        @return: List of feature IDs
        """
        result = []
        if ftype not in self.feat2id:
            self.feat2id[ftype] = {}
        fdict = self.feat2id[ftype]
        for feat in feats:
            if feat not in fdict:
                featid = len(self.feats)
                self.feats.append((feat,ftype))
                fdict[feat] = featid
            result.append(fdict[feat])
        return result

class TermCounts(dict):
    """Mapping from feature ID's to number of occurrences.

    @ivar docs: Number of documents processed
    @ivar total: Total occurrences (== sum of values)
    @ivar filename: File to store feature counts (may be None)
    """

    def __new__(cls, filename=None, items=None):
        """For loading the counts from a pickle instead"""
        if isinstance(filename,path) and filename.exists():
            return cPickle.load(file(filename, "rb"))
        else:
            return super(TermCounts, cls).__new__(cls, filename)

    def __init__(self, filename=None, items=None):
        """Initialise the term counter to zero.

        @param filename: Optional path object to a pickle with the
        instance.
        @param items: Optional list of items, each of which is given
        to the add method.
        """
        if not (isinstance(filename,path) and filename.exists()):
            self.docs = 0
            self.total = 0
        if items is not None:
            for item in items:
                self.add(item)
        self.filename = filename

    def dump(self):
        """Write a pickle of a TermCounts instance to disk, keeping
        temporary backup"""
        if self.filename is None:
            return
        makeBackup(self.filename)
        cPickle.dump(self, file(self.filename, "wb"), protocol=2)
        removeBackup(self.filename)

    def add(self, features):
        """Add a list of features to cumulative feature frequencies"""
        self.docs += 1
        self.total += len(features)
        for f in features:
            if f not in self:
                self[f] = 1
            else:
                self[f] += 1

    def subtract(self, other):
        """Return a TermCounts instance representing the article from
        this one less the articles from other."""
        result = TermCounts()
        result.docs = self.docs - other.docs
        result.total = self.total - other.total
        for termid,count in self.iteritems():
            result[termid] = count - other.get(termid,0)
        return result
        
