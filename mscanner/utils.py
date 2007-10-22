"""Dumping ground for utility functions that don't naturally fit anywhere and
are used in several modules."""

from __future__ import with_statement
from __future__ import division


                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


def update(obj, vars, exclude=['self']):
    """Update instance attributes (using a dictionary)
    
    Example: C{update(self, locals())}
    
    @param obj: Instance to update using setattr()
    @param vars: Dictionary of variables to store
    @param exclude: Variables to exclude, defaults to ['self']
    """
    for k, v in vars.iteritems():
        if k not in exclude:
            setattr(obj, k, v)


def delattrs(obj, *vars):
    """Remove named attributes from instance
    
    Example: C{delattrs(self, "_property", "_other")}
    
    @param objs: Object to update via delattr
    @param vars: Instance attribute names to delete
    """
    for ivar in vars:
        try:
            delattr(obj, ivar)
        except AttributeError:
            pass


def make_random_subset(k, pool, exclude):
    """Choose a random subset of k articles from pool

    Use this function when the pool is large (say, 16 million items), we don't
    mind if the order of pool gets scrambled in the process, and we need to
    exclude certain items from being selected.
    
    @param k: Number of items to choose from pool
    @param pool: Array of items to choose from (will be scrambled!)
    @param exclude: Set of items that may not be chosen
    @return: A new array of chosen items
    """
    from random import randint
    import numpy as nx
    n = len(pool)
    assert 0 <= k <= n
    for i in xrange(k):
        # Non-selected items are in 0 ... n-i-1
        # Selected items are n-i ... n
        dest = n-i-1
        choice = randint(0, dest) # 0 ... n-i-1 inclusive
        while pool[choice] in exclude:
            choice = randint(0, dest)
        # Move the chosen item to the end, where so it will be part of the
        # selected items in the next iteration. Note: this works using single
        # items - it but would break with slices due to their being views into
        # the vector.
        pool[dest], pool[choice] = pool[choice], pool[dest]
    # Phantom iteration: selected are n-k ... n
    return nx.array(pool[n-k:])


def usetempfile(function):
    """Decorator to call a method with a temporary file
    
    Create a temporary file, pass the file path to the wrapped function and
    finally remove the file afterwards. Meant for wrapping unit testing methods
    which require access to a temporary file."""
    import tempfile
    from path import path
    def tempfile_wrapper(self):
        try:
            fpath = path(tempfile.mktemp())
            return function(self, fpath)
        finally:
            if fpath.isfile():
                fpath.remove()
    return tempfile_wrapper


class FileTransaction(file):
    """Transaction for Cheetah templates to output direct-to-file.
    
    Cheetah defaults to DummyTransaction which creates a huge list and
    joins them up to create a string.  This is way slower than writing to
    file directly.
    
    Usage::
        with FileTransaction("something.html","wb") as ft:
            Template().respond(ft)
    """
    
    def __init__(self, *args, **kw):
        """Open the file, same parameters as for the builtin"""
        file.__init__(self, *args, **kw)
        self.response = self

    def writeln(self):
        """Write a line of output"""
        self.write(txt)
        self.write('\n')

    def getvalue(self):
        """Not implemented"""
        return None

    def __call__(self):
        return self


def read_pmids(filename, 
              include=None, 
              exclude=None, 
              broken_name=None, 
              exclude_name=None, 
              withscores=False):
    """Read PubMed IDs one per line from filename.

    @param filename: Path to file with PubMed IDs, formatted one PubMed ID per
    line, with optional score after the PubMed ID. Blank lines and lines
    starting with # are ignored, as is data after the PMID and score.
    
    @param include: Only return members of this set (other PubMed IDs are
    considered "broken").

    @param broken_name: File to write non-included ("broken") PubMed IDs
    
    @param exclude: Do not return members of this set
    
    @param exclude_name: File to write excluded PubMed IDs
    
    @param withscores: Also read the score after the PubMed ID on each line.
    
    @returns: Iterator over PubMed ID, or (Score, PubMed ID) if withscores
    True. """
    import logging as log
    count = 0
    broken = []
    excluded = []
    if filename.isfile(): 
        with open(filename, "r") as infile:
            for line in infile:
                sline = line.strip()
                if sline == "" or sline.startswith("#"):
                    continue
                splits = sline.split()
                pmid = int(splits[0])
                if include is not None and pmid not in include:
                    broken.append(pmid)
                    continue
                if exclude is not None and pmid in exclude:
                    excluded.append(pmid)
                    continue
                if withscores:
                    yield float(splits[1]), pmid
                else:
                    yield pmid
                count += 1
    if broken_name is not None:
        with open(broken_name, "w") as f:
            f.write("\n".join(str(s) for s in broken))
    if exclude_name is not None:
        with open(exclude_name, "w") as f:
            f.write("\n".join(str(s) for s in excluded))
    if count == 0:
        raise ValueError("No valid PMIDs found in file")
    log.debug("Got %d PubMed IDs from %s", count, filename.basename())

    
def write_scores(filename, pairs):
    """Write scores and PubMed IDs to file, in decreasing order of score.
    @param pairs: Iterable of (score, PMID)     
    """
    from path import path
    path(filename).write_lines(
        "%-10d %f" % (p,s) for s,p in sorted(pairs, reverse=True))
    
    
def load_articles(article_db_path, pmidlist_path):
    """Retrieve Article objects given a file of PubMed IDs.

    @param article_db_path: Path to a berkeley DB mapping PubMed IDs
    to Article objects.

    @param pmidlist_path: Path to a text file listing one PubMed ID per line.

    @return: List of Article objects in the order given in the text file.

    @note: The first called with a given PMID list caches the results in a
    .pickle, and later calls load the pickle."""
    import cPickle
    from mscanner.medline import Shelf
    from contextlib import closing
    from path import path # used in the line below
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        with open(cache_path, "rb") as f:
            return cPickle.load(f)
    pmids = read_pmids(pmidlist_path)
    with closing(Shelf.open(article_db_path, "r")) as artdb:
        articles = [artdb[str(p)] for p in pmids]
    with open(cache_path, "wb") as f:
        cPickle.dump(articles, f, protocol=2)
    return articles


def no_valid_pmids_page(filename, pmids):
    """Print an error page when no valid PMIDs were found
    @param filename: Path to output file
    @param pmids: List of any provided PMIDs (all invalid)
    """
    from Cheetah.Template import Template
    from mscanner import utils
    from mscanner.configuration import rc
    import logging as log
    log.warning("No valid PubMed IDs were found!")
    with utils.FileTransaction(filename, "w") as ft:
        page = Template(file=str(rc.templates/"notfound.tmpl"))
        page.notfound_pmids = pmids
        page.respond(ft)
