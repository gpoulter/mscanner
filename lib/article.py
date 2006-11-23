"""Handle Articles, processed files, features, feature occurrence counts

@author: Graham Poulter
                                     

createStatusFile() -- Start a status file
createStatusFile() -- Update progress in a status file
chooseRandomLines() -- Get random lines from a file
readPMIDFile() -- Get PMIDs reliably from a text file
getArticles() -- Retrieve Articles from a DB, caching results in a Pickle
makeBackup() -- Create a .recovery file
removeBackup() -- Remove a .recovery file

Article -- Stores attributes of an article
FileTracker -- Track processed files (to avoid re-parsing), with on-disk persistence
FeatureMapping -- Mapping between string features and 16-bit feature IDs
TermCounts -- Store number of occurrences of each feature, + total occurrences and documents

"""

from array import array
import cPickle
import dbshelve
import logging as log
import os
from path import path
import time

def createStatusFile(statfile, dataset, total=0):
    """Create a file containing PID, 0, 0, dataset, str(time.time())"""
    if statfile.exists():
        raise RuntimeError("MedScanner instance is already running")
    statfile.write_text("%d\n%d\n%d\n%s\n%s\n" % (os.getpid(), 0, 0, dataset, str(time.time())))

def updateStatusFile(statfile, progress, total=None):
    """Update the 0, 0 lines with progress and (optional) total.  If
    progress is None, set to total."""
    lines = statfile.lines()
    if total is not None:
        lines[2] = str(total)
    if progress is not None:
        lines[1] = str(progress)
    else:
        lines[1] = lines[2]
    statfile.write_lines(lines)

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
        if meshterms is None: self.meshterms = set()
        self.authors = authors
        if authors is None: self.authors = set()
        self.chemicals = chemicals
        if chemicals is None: self.chemicals = set()

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
        self.trackfile = path(trackfile)
        self.load()

    def dump(self):
        """Write the list of tracked files, one per line"""
        if self.trackfile == "None": return
        proclist = list(self)
        proclist.sort()
        makeBackup(self.trackfile)
        file(self.trackfile, "w").write("\n".join(proclist))
        removeBackup(self.trackfile)

    def load(self):
        """Read the list of files into the tracker"""
        self.clear()
        if self.trackfile.isfile():
            self.update(line.strip() for line in self.trackfile.lines())

    def add(self, fname):
        """Add a file to the tracker"""
        set.add(self, fname.basename())

    def toprocess(self, files):
        """Given a list of files return sorted list of those not in the tracker"""
        toprocess = [ f for f in files if f.basename() not in self ]
        toprocess.sort()
        return toprocess

class FeatureMapping:
    """Curates a database of encountered term features, providing
    methods to map between terms and 16-bit integer IDs.

    @ivar termid: A dict mapping string to ID
    @ivar term: A list mapping ID to string
    """

    def __init__(self, termfile = None):
        """Initialise the database
        @param termfile: Path to text file with list of terms
        """
        self.termfile = path(termfile)
        self.termid = {}
        self.term = []
        self.load()

    def load(self):
        """Load the term mapping from disk"""
        self.termid = {}
        self.term = []
        if not self.termfile.exists(): return
        f = file(self.termfile, "r")
        for term in f:
            term = term.strip()
            self.termid[term] = len(self.term)
            self.term.append(term)
        f.close()

    def dump(self):
        """Write the term mapping to disk"""
        if self.termfile == "None": return
        makeBackup(self.termfile)
        f = file(self.termfile, "w")
        for term in self.term:
            f.write(term+"\n")
        f.close()
        removeBackup(self.termfile)

    def __getitem__(self, feature_id):
        """Return feature string given feature ID"""
        return self.term[feature_id]

    def __len__(self):
        """Return number of distinct features"""
        return len(self.term)

    def getterms(self, term_ids):
        """Return feature strings given feature ids"""
        return [ self.term[tid] for tid in term_ids ]

    def getids(self, terms):
        """Get term IDs corresponding to the given list of terms.
        @return: array('H') of term IDs
        @note: Dynamically creates new term IDs as necessary.
        """
        result = []
        for term in terms:
            if term not in self.termid:
                tid = len(self.term)
                self.term.append(term)
                self.termid[term] = tid
            result.append(self.termid[term])
        return array("H",result)

class TermCounts(dict):
    """Mapping from feature ID's to number of occurrences.

    @ivar docs: Number of documents processed
    @ivar total: Total occurrences (== sum of values)
    """
    __slots__ = [ "docs", "total" ]

    def __init__(self, items=None):
        """Initialise the term counter to zero.
        @note: If items is a list of features vectors, add those
        feature vectors to this instance.
        """
        dict.__init__(self)
        self.docs = 0
        self.total = 0
        if items is not None:
            for features in items:
                self.add(features)         

    @staticmethod
    def load(filepath):
        """Return TermCounts instance, either loaded from pickle or
        freshly instantiated""" 
        if path(filepath).exists():
            return cPickle.load(file(filepath, "rb"))
        else:
            return TermCounts()

    @staticmethod
    def dump(instance, filepath):
        """Write a pickle of a TermCounts instance to disk, keeping
        temporary backup"""
        filepath = path(filepath)
        if filepath == "None": return
        makeBackup(filepath)
        cPickle.dump(instance, file(filepath, "wb"), protocol=2)
        removeBackup(filepath)

    def add(self, features):
        """Adds the terms from an article to the term counts"""
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
