#!/usr/bin/env python

"""Perform queries with a gene-drug co-occurrence filtering step

Choose operation by executing a Python expression::

    python pharmdemo/main.py 'pharmdemo("pharmgkb-070205.txt")'
    
Which will perform a query using PharmGKB PubMed IDs, then filter
the results for gene-drug occurrences, and export a database for
used on http://pharmdemo.stanford.edu
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

from path import path
import logging as log
import sys

from mscanner.configuration import rc, initLogger
from mscanner import scorefile
from mscanner import queryenv
from mscanner import validenv
from pharmdemo import genedrug
from pharmdemo import dbexport

class PharmdemoQuery(queryenv.QueryEnvironment):

    def getGDFilterResults(self, pmids_path, export_db=False):
        """Filter L{results} and L{inputs} for those gene-drug co-occurrences
        
        @param pmids_path: Path to input PubMed IDs
        
        @param export_db: If True, export associations for PharmDemo
        
        @return: Set of articles with gene-drug associations.
        """
        from itertools import chain
        self.input_pmids = set(scorefile.readPMIDs(pmids_path, include=self.featdb))
        self.featinfo = self.getFeatureInfo()
        self.getResults()
        log.debug("Gene-drug associations on results")
        gdfilter = genedrug.getGeneDrugFilter(
            rc.genedrug_db, rc.drugtable, rc.gapscore_db)
        gdarticles = []
        for score, pmid in chain(self.results, self.inputs):
            a = self.artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                gdarticles.append(a)
        gdfilter.close()
        if export_db == True:
            log.debug("Exporting database")
            gdexport = dbexport.GeneDrugExport(gdarticles)
            gdexport.writeGeneDrugCountsCSV(rc.genedrug_csv)
            gdexport.exportText(rc.genedrug_sql)
        return gdarticles



class PharmdemoValidation(validenv.ValidationEnvironment):
    
    def genedrug_filter(self):
        """Create a membership test for gene-drug association.
        
        To use the filter, assign this method to self.postFilterFunction
        
        @return: Set of PubMed IDs which have gene-drug co-occurrences.
        """
        cwd = rc.valid_report_dir
        from pharmdemo import genedrug 
        log.info("Getting gene-drug associations") 
        pos_arts = scorefile.getArticles(cwd/rc.articledb, self.positives)
        neg_arts = scorefile.getArticles(cwd/rc.articledb, self.negatives)
        gdfilter = genedrug.getGeneDrugFilter(
            cwd/rc.genedrug, cwd/rc.drugtable, cwd/rc.gapscore)
        postfilter = set()
        for art in chain(pos_arts, neg_arts):
            gdresult = gdfilter(art)
            art.genedrug = gdresult
            if len(gdresult) > 0:
                postfilter.add(art.pmid)
        return postfilter



def pharmdemo(pmidfile):
    """Perform a standard query, then filter the Query which exports to the PharmDemo database for PharmGKB
    
    @param pmidfile: Name of the file under rc.corpora to use for input PubMed IDs
    """
    rc.dataset = "pharmdemo"
    rc.threshold = 20.0 # Higher than usual threshold
    rc.limit = 10000 # Want lots of results
    env = PharmdemoQuery()
    env.getGDFilterResults(rc.corpora / pmidfile, export_db=True)


if __name__ == "__main__":
    initLogger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression"
    else:
        eval(sys.argv[1])
