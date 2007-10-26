"""I/O functions - for reading and writing certain file formats."""

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

    
def write_scores(filename, pairs, sort=True):
    """Write scores and PubMed IDs to file, in decreasing order of score.
    @param pairs: Iterable of (score, PMID)     
    """
    from path import path
    spairs = sorted(pairs, reverse=True) if sort else pairs
    path(filename).write_lines(
        "%-10d %f" % (p,s) for s,p in spairs)


def no_valid_pmids_page(filename, pmids):
    """Print an error page when no valid PMIDs were found
    @param filename: Path to output file
    @param pmids: List of any provided PMIDs (all invalid)
    """
    from Cheetah.Template import Template
    from mscanner.configuration import rc
    import logging as log
    log.warning("No valid PubMed IDs were found!")
    with FileTransaction(filename, "w") as ft:
        page = Template(file=str(rc.templates/"notfound.tmpl"))
        page.notfound_pmids = pmids
        page.respond(ft)
