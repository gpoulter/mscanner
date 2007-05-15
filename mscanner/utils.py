"""Handle Articles, processed files, features, feature occurrence counts

preserve_cwd -- Decorator to preserve working directory
countFeatures() -- Count occurrences of features in a set of articles
readPMIDs() -- Get PMIDs (and scores too) from a text file
writePMIDScores() -- Write PMIDs and scores to a text file
runMailer() -- Send e-mails to people listed in a file
getArticles() -- Retrieve Articles from a DB, caching results in a Pickle
FileTracker -- Track processed files (to avoid re-parsing), with on-disk persistence

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import cPickle
import logging as log
import numpy as nx
import os
from path import path

from mscanner import dbshelve

def preserve_cwd(f):
    """Decorator which saves working directory before callijng the function,
    and restores it afterwards."""
    def cwd_preserver(*a, **kw):
        try:
            cwd = os.getcwd()
            return f(*a, **kw)
        finally:
            os.chdir(cwd)
    def decorator_update(dest, src):
        dest.__name__ = src.__name__
        dest.__doc__ = src.__doc__
        dest.__dict__.update(src.__dict__)
    decorator_update(cwd_preserver, f)
    return cwd_preserver

def selfupdate(onlyargs=False, exclude=[]):
    """Call in any method to set instance attributes from local variables.
    
    @param onlyargs: If True, use only named arguments
    
    @param exclude: Names of other variables to exclude.

    @note: Instance to update is assumed to be first argument of the caller.
    
    @note: Equivalent to vars(self).update(vars()); del self.self
    """
    import inspect
    # Get caller frame (must be disposed of!)
    cf = inspect.currentframe().f_back 
    try:
        # Instance is first argument to the caller
        instance = cf.f_locals[cf.f_code.co_varnames[0]]
        # Get names of variables to use
        if onlyargs:
            varnames = cf.f_code.co_varnames[1:cf.f_code.co_argcount]
        else:
            varnames = cf.f_code.co_varnames[1:]
        # Update instance with those names
        for var in varnames:
            if var not in exclude:
                setattr(instance, var, cf.f_locals[var])
    finally:
        del cf
        
def update(instance, variables, exclude=['self']):
    """Update instance attributes
    
    @param instance: Instance to update via setattr()
    @param variables: Dictionary of variables
    @param exclude: Variables to exclude, defaults to ['self']
    """
    for k, v in variables.iteritems():
        if k not in exclude:
            setattr(instance, k, v)

def readPMIDs(filename, include=None, exclude=None, withscores=False):
    """Read PubMed IDs one per line from filename.

    @param filename: Path to file containing one PubMed ID per line, with
    optional score after the PubMed ID. File format ignores blank lines and
    lines starting with #, and only parses the line up to the first whitespace
    character.

    @param include: Only return members of this set

    @param exclude: Do not return members of this set
    
    @param withscores: Also read the score after the PubMed ID on each line.
    
    @returns: An iterator over PubMed ID, or 
    (Score, PubMed ID) if withscores==True
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
            yield float(splits[1]), pmid
        else:
            yield pmid
        count += 1
    if count == 0:
        raise RuntimeError("Did not succeed in reading any PMIDs from %s" %
                           filename)

def writePMIDScores(filename, pairs):
    """Write (score, PMID) pairs to file, in decreasing order of score."""
    filename.write_lines("%-10d %f" % (p,s) for s,p in \
                         sorted(pairs, reverse=True))

def runMailer(smtpserver, mailer):
    """Read e-mail addresses from mail file and send them a message saying
    MScanner is available.

    @param mailer: Path object to file of e-mail addresses
    """
    if not mailer.isfile():
        return
    import smtplib
    import logging
    server = smtplib.SMTP(smtpserver)
    server.set_debuglevel(0)
    fromaddr = "nobody@mscanner.stanford.edu"
    for email in mailer.lines(retain=False):
        logging.debug("Sending availability alert to %s", email)
        msg = "From: %s\r\nTo: %s\r\n\r\n" % (fromaddr, email)
        msg += """
Someone requested that an alert be sent to this address when MScanner,
a service for classifying Medline citations given a set of examples
(http://mscanner.stanford.edu) completed its current request.

MScanner has completed the last submission and at the time of sending
this email is available for query or validation operations.

This is a once-off notification. The list of e-mails to notify
is erased immediately after sending.
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

    @note: The first time it is called for a given pmidlist_path, the results
    are cached in a .pickle appended to pmidlist_path, and later calls simply
    use the cached results. 
    """
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        return cPickle.load(file(cache_path, "rb"))
    pmids = readPMIDs(pmidlist_path)
    artdb = dbshelve.open(article_db_path, 'r')
    articles = [artdb[str(p)] for p in pmids]
    cPickle.dump(articles, file(cache_path,"wb"), protocol=2)
    artdb.close()
    return articles

class FileTracker(set):
    """Class which tracks processed files.

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

