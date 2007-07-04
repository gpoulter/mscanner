"""Handle reading and writing of PMIDs and their scores

readPMIDs() -- Get PMIDs (and scores too) from a text file
writePMIDScores() -- Write PMIDs and scores to a text file
getArticles() -- Retrieve Articles from a DB, caching results in a Pickle
readDescriptor() -- Read paramaters from a descriptor file
writeDescriptor() -- Write paramaters to a descriptor file

                                   
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

def readPMIDs(filename, outbase=None, include=None, exclude=None, withscores=False):
    """Read PubMed IDs one per line from filename.

    @param filename: Path to file containing one PubMed ID per line, with
    optional score after the PubMed ID. File format ignores blank lines and
    lines starting with #, and only parses the line up to the first whitespace
    character.

    @param include: Only return members of this set

    @param exclude: Do not return members of this set
    
    @param withscores: Also read the score after the PubMed ID on each line.
    
    @returns: An iterator over PubMed ID, or (Score, PubMed ID) if
    withscores is True.
    """
    import logging as log
    count = 0
    broken_file = file(outbase+".broken.txt","w") if outbase else None
    exclude_file = file(outbase+".exclude.txt","w") if outbase else None
    for line in file(filename, "r"):
        sline = line.strip()
        if sline == "" or sline.startswith("#"):
            continue
        splits = sline.split()
        pmid = int(splits[0])
        if include is not None and pmid not in include:
            if outbase: broken_file.write(str(pmid)+"\n")
            log.warn("PMID %d is not a member of include" % pmid)
            continue
        if exclude is not None and pmid in exclude:
            if outbase: exclude_file.write(str(pmid)+"\n")
            log.warn("PMID %d is a member of exclude" % pmid)
            continue
        if withscores:
            yield float(splits[1]), pmid
        else:
            yield pmid
        count += 1
    if outbase:
        broken_file.close()
        exclude_file.close()
    if count == 0:
        raise ValueError("Did not succeed in reading any non-excluded PMIDs"+\
                         "from %s" % filename)
    else:
        log.debug("Got %d PubMed IDs from %s", count, filename.basename())

def writePMIDScores(filename, pairs):
    """Write (score, PMID) pairs to filename, in decreasing
    order of score."""
    from path import path
    path(filename).write_lines(
        "%-10d %f" % (p,s) for s,p in sorted(pairs, reverse=True))
    
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
    import cPickle
    from path import path
    from mscanner.support import dbshelve
    cache_path = path(pmidlist_path + ".pickle")
    if cache_path.isfile():
        return cPickle.load(file(cache_path, "rb"))
    pmids = readPMIDs(pmidlist_path)
    artdb = dbshelve.open(article_db_path, 'r')
    articles = [artdb[str(p)] for p in pmids]
    cPickle.dump(articles, file(cache_path,"wb"), protocol=2)
    artdb.close()
    return articles

descriptor_keys = dict(
    alpha=float,
    code=str,
    dataset=str,
    limit=int,
    nfolds=int,
    numnegs=int,
    operation=str,
    threshold=float,
    timestamp=float,
)
    
def readDescriptor(fpath):
    """Reads a descriptor file, returning a dictionary of parameters.

    @note: Each line is formatted as "#key = value".  The reading stops
    at the first line that does not start with '#'.   The valid keys
    are listed in the converter dictionary.

    @note: The same descriptor file can be used as a PubMed ID list,
    as the PubMed-ID reader ignores lines beginning with '#'.
    """

    f = file(fpath, "r")
    line = f.readline()
    from mscanner.support.storage import Storage
    result = Storage()
    while line.startswith("#"):
        key, value = line[1:].split(" = ",1)
        value = descriptor_keys[key](value.strip())
        result[key] = value
        line = f.readline()
    f.close()
    return result
    
def writeDescriptor(fpath, pmids, params):
    """Write parameters and PubMed IDs to the descriptor file.
    
    @param fpath: File to write
    @param pmids: List of PubMed IDs, may be None
    @param params: Key-value dictionary to write. Values are converted with
    str().  Only certain keys are written."""
    f = file(fpath, "w")
    for key, value in params.iteritems():
        if key not in descriptor_keys: 
            continue
        f.write("#" + key + " = " + str(value).strip() + "\n")
    if pmids is not None:
        fpath.write_lines([str(pmid) for pmid in pmids])
    f.close()
