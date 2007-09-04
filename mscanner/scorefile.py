"""Performs reading/writing of various file types

Functions in this file:
    - L{readPMIDs}: Get PMIDs (and scores too) from a text file
    - L{writePMIDScores}: Write PMIDs and scores to a text file
    - L{getArticles}: Retrieve Articles from a DB, caching results in a Pickle
    - L{readDescriptor}: Read paramaters from a descriptor file
    - L{writeDescriptor}: Write paramaters to a descriptor file
"""

                                     
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


def readPMIDs(filename, include=None, exclude=None, 
              broken_name=None, exclude_name=None, withscores=False):
    """Read PubMed IDs one per line from filename.

    @param filename: Path to file with PubMed IDs. Format is one PubMed ID per
    line, with optional score after the PubMed ID. File format ignores blank
    lines and lines starting with #, and only parses the line up to the first
    whitespace character.
    
    @param include: Only return members of this set (other PubMed IDs are
    considered "broken").

    @param broken_name: File to write non-included ("broken") PubMed IDs
    
    @param exclude: Do not return members of this set
    
    @param exclude_name: File to write excluded PubMed IDs
    
    @param withscores: Also read the score after the PubMed ID on each line.
    
    @returns: Iterator over PubMed ID, or (Score, PubMed ID) if withscores is
    True. """
    import logging as log
    count = 0
    broken = []
    excluded = []
    for line in open(filename, "r"):
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
    if broken_name != None:
        broken_file = open(broken_name, "w")
        broken_file.write("\n".join(str(s) for s in broken))
        broken_file.close()
    if exclude_name != None:
        exclude_file = open(exclude_name, "w")
        exclude_file.write("\n".join(str(s) for s in excluded))
        exclude_file.close()
    log.debug("Got %d PubMed IDs from %s", count, filename.basename())


def writePMIDScores(filename, pairs):
    """Write scores and PubMed IDs to file, in decreasing order of score.
    
    @param pairs: Iterable of (score, PMID)     
    """
    from path import path
    path(filename).write_lines(
        "%-10d %f" % (p,s) for s,p in sorted(pairs, reverse=True))
    
    
def getArticles(article_db_path, pmidlist_path):
    """Retrieve Article objects given a file of PubMed IDs.

    @param article_db_path: Path to a berkeley DB mapping PubMed IDs
    to Article objects.

    @param pmidlist_path: Path to a text file listing one PubMed ID per line.

    @return: List of Article objects in the order given in the text file.

    @note: The first called with a given PMID list caches the results in a
    .pickle, and later calls load the pickle."""
    import cPickle
    from path import path
    from mscanner.support import dbshelve
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        return cPickle.load(open(cache_path, "rb"))
    pmids = readPMIDs(pmidlist_path)
    artdb = dbshelve.open(article_db_path, 'r')
    articles = [artdb[str(p)] for p in pmids]
    cPickle.dump(articles, open(cache_path,"wb"), protocol=2)
    artdb.close()
    return articles


def parsebool(s):
    """Handler for converting strings to booleans"""
    if isinstance(s, basestring):
        s = s.strip()
        if s == "0" or s == "False":
            return False
        elif s == "1" or s == "True":
            return True
        else:
            raise ValueError("Failed to parse boolean: %s" % s)
    else:
        return bool(s)


descriptor_keys = dict(
    alpha=float,
    code=str,
    dataset=str,
    delcode=str,
    hidden=parsebool,
    limit=int,
    nfolds=int,
    numnegs=int,
    operation=str,
    threshold=float,
    timestamp=float,)


def readDescriptor(fpath):
    """Reads a descriptor file, returning a dictionary of parameters.

    Each descriptor line is formatted as "#key = value". Reading stops at the
    first line that does not start with '#'. Valid keys are in
    L{descriptor_keys}. The same file can be used with L{readPMIDs}, which
    ignores lines beginning with '#'.

    @return: Storage object, with additional "_filename" key that contains
    fpath."""
    f = open(fpath, "r")
    line = f.readline()
    from mscanner.support.storage import Storage
    result = Storage()
    while line.startswith("#"):
        key, value = line[1:].split(" = ",1)
        value = descriptor_keys[key](value.strip())
        result[key] = value
        line = f.readline()
    result["_filename"] = fpath
    f.close()
    return result


def writeDescriptor(fpath, pmids, params):
    """Write parameters and PubMed IDs to the descriptor file.
    
    @param fpath: File to write

    @param pmids: List of PubMed IDs, may be None
    
    @param params: Key-value dictionary to write. Values are converted with
    str(). Only keys that have a member in descriptor_keys are actually written
    to the file."""
    f = open(fpath, "w")
    for key, value in params.iteritems():
        if key in descriptor_keys: 
            f.write("#" + key + " = " + str(value) + "\n")
    if pmids is not None:
        for pmid in pmids:
            f.write(str(pmid)+"\n")
    f.close()


def no_valid_pmids_page(pmids):
    """Print a page when no valid PMIDs were found
    
    The page includes a list of any submitted IDs and links to PubMed (
    however, these IDs were not part of MScanner's database). """
    import logging
    from mscanner.support.gcheetah import TemplateMapper, FileTransaction
    from mscanner.configuration import rc
    logging.warning("No valid PubMed IDs were found!")
    mapper = TemplateMapper(root=rc.templates, kwargs=dict(filter="Filter"))
    ft = FileTransaction(rc.report_index, "w")
    index = mapper.notfound(
        rc = rc,
        notfound_pmids = pmids).respond(ft)
    ft.close()

