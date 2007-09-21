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

from itertools import chain
import logging as log
import sys

from mscanner.configuration import rc, initLogger
from mscanner import scorefile, queryenv, validenv
from pharmdemo import genedrug, dbexport


### GENE-DRUG CONFIGURATION
rc.genedrug = lambda: rc.working / "genedrug"
## Path to GAPScore results db
rc.gapscore_db = lambda: rc.genedrug / "gapscore.db"
## Path to gene-drug co-occurrence cache
rc.genedrug_db = lambda: rc.genedrug / "genedrug.db"
## Path to drug table
rc.drugtable = lambda: rc.genedrug / "drugtable.txt"
## Path to CSV table of co-occurrences
rc.genedrug_csv = lambda: rc.genedrug / "pharmdemo.csv"
## Path to SQL table of co-occurrences
rc.genedrug_sql = lambda: rc.genedrug / "pharmdemo.sql"


class PharmdemoQuery(queryenv.Query):

    def genedrug_query(self, input, export_db=False):
        """Filter L{results} and L{inputs} for those gene-drug co-occurrences
        
        @param filename: Path to input PubMed IDs
        
        @param export_db: If True, export associations for PharmDemo
        
        @return: Set of articles with gene-drug associations.
        """
        import time
        if rc.timestamp is None: 
            rc.timestamp = time.time() 
        self.load_pmids(input)
        if len(self.pmids) == 0: return
        self.make_featinfo()
        self.make_results()
        # Do associations
        log.debug("Gene-drug associations on results")
        gdfilter = genedrug.getGeneDrugFilter(
            rc.genedrug_db, rc.drugtable, rc.gapscore_db)
        gd_articles = []
        gd_pmids = set()
        for score, pmid in chain(self.results, self.inputs):
            a = self.env.artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                gd_articles.append(a)
                gd_pmids.add(pmid)
        gdfilter.close()
        self.inputs  = [ (s,p) for s,p in self.inputs  if p in gd_pmids ]
        self.results = [ (s,p) for s,p in self.results if p in gd_pmids ]
        if export_db == True:
            log.debug("Exporting database")
            gdexport = dbexport.GeneDrugExport(gd_articles)
            gdexport.writeGeneDrugCountsCSV(rc.genedrug_csv)
            gdexport.exportText(rc.genedrug_sql)
        # Finish by writing results
        self.save_results()
        self.write_report()
        rc.timestamp = None # reset for next run
        log.info("FINISHING QUERY %s", rc.dataset)



class PharmdemoValidation(validenv.Validation):
    
    def genedrug_filter(self):
        """Create a membership test for gene-drug association.
        
        To use the filter, assign this method to self.postFilterFunction
        
        @return: Set of PubMed IDs which have gene-drug co-occurrences.
        """
        log.info("Getting gene-drug associations") 
        pos_arts = scorefile.load_articles(rc.articledb, self.positives)
        neg_arts = scorefile.load_articles(rc.articledb, self.negatives)
        gdfilter = genedrug.getGeneDrugFilter(
            rc.genedrug, rc.drugtable, rc.gapscore)
        postfilter = set()
        for art in chain(pos_arts, neg_arts):
            gdresult = gdfilter(art)
            art.genedrug = gdresult
            if len(gdresult) > 0:
                postfilter.add(art.pmid)
        return postfilter



def pharmdemo():
    """Perform a standard query, then filter the Query which exports to the PharmDemo database for PharmGKB
    
    @param pmidfile: Name of the file under rc.corpora to use for input PubMed IDs
    """
    filename = rc.corpora / "pharmgkb-070205.txt"
    #filename = rc.corpora / "genedrug-small.txt"
    rc.dataset = "pharmdemo"
    rc.threshold = 20.0 # Higher than usual threshold
    rc.limit = 10000 # Want lots of results
    op = PharmdemoQuery(rc.working / "query" / rc.dataset)
    op.genedrug_query(filename, export_db=True)


if __name__ == "__main__":
    initLogger()
    pharmdemo()
