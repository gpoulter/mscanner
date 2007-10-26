#!/usr/bin/env python

"""Perform queries with a gene-drug co-occurrence filtering step

Choose operation by executing a Python expression::

    python pharmdemo/main.py 'pharmdemo("pharmgkb-070205.txt")'
    
Which will perform a query using PharmGKB PubMed IDs, then filter
the results for gene-drug occurrences, and export a database for
used on http://pharmdemo.stanford.edu
"""


from itertools import chain
import logging as log
import sys

from mscanner.configuration import rc, start_logger
from mscanner.core.QueryManager import QueryManager
from mscanner.core.ValidationManager import ValidationManager
from mscanner.medline import Databases
from mscanner.pharmdemo import genedrug
from mscanner.pharmdemo.Exporter import Exporter


                                     
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


class PharmdemoQuery(QueryManager):

    def genedrug_query(self, input, export_db=False):
        """Query where L{results} and L{inputs} are filtered for gene-drug
        co-occurrences, and the associations are exported as an SQL file.
        
        @param input: Path to PubMed IDs
        
        @param export_db: If True, export associations for PharmDemo
        """
        if not self._load_input(input):
            return
        self._make_feature_info()
        self._make_results()
        # Carry out associations
        log.debug("Gene-drug associations on results")
        gdfinder = genedrug.open_genedrug_finder(
            rc.genedrug_db, rc.drugtable, rc.gapscore_db)
        gd_articles = []
        gd_pmids = set()
        for score, pmid in chain(self.results, self.inputs):
            a = self.env.artdb[str(pmid)]
            a.genedrug = gdfinder[a]
            if len(a.genedrug) > 0:
                gd_articles.append(a)
                gd_pmids.add(pmid)
        gdfinder.close()
        self.inputs  = [ (s,p) for s,p in self.inputs  if p in gd_pmids ]
        self.results = [ (s,p) for s,p in self.results if p in gd_pmids ]
        if export_db == True:
            log.debug("Exporting database")
            gdexport = Exporter(gd_articles)
            gdexport.write_genedrug_csv(rc.genedrug_csv)
            gdexport.export_sqlfile(rc.genedrug_sql)
        # Finish by writing results
        self._save_results()
        self._write_report()



class PharmdemoValidation(ValidationManager):
    
    def genedrug_filter(self):
        """Membership test for gene-drug association.
        
        @return: Set of PubMed IDs which have gene-drug co-occurrences.
        """
        log.info("Getting gene-drug associations") 
        pos_arts = Databases.load_articles(rc.articledb, self.positives)
        neg_arts = Databases.load_articles(rc.articledb, self.negatives)
        gdfinder = genedrug.open_genedrug_finder(
            rc.genedrug, rc.drugtable, rc.gapscore)
        postfilter = set()
        for art in chain(pos_arts, neg_arts):
            gdresult = gdfinder[art]
            art.genedrug = gdresult
            if len(gdresult) > 0:
                postfilter.add(art.pmid)
        return postfilter



def pharmdemo():
    """Perform a query with exporting of Pharmdemo associations"""
    #filename = rc.corpora / "pharmgkb-070205.txt"
    filename = rc.corpora / "genedrug-small.txt"
    rc.dataset = "pharmdemo"
    rc.threshold = 20.0 # Higher than usual threshold
    rc.limit = 10000 # Want lots of results
    op = PharmdemoQuery(rc.working / "query" / rc.dataset)
    op.genedrug_query(filename, export_db=True)


if __name__ == "__main__":
    start_logger()
    pharmdemo()
