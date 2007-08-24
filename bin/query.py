#!/usr/bin/env python

"""Query the MEDLINE database with a set of positive articles

Usage as a backend to the web interface:

query.py dataset limit threshold positives_path

For experiments, use a Python snippet:

query.py 'query("pg04","pg07")'

                                   
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

from path import path
import sys

from mscanner.configuration import rc, initLogger
from mscanner.queryenv import QueryEnvironment
            
ds_map = {
    "aids"        : "aids-bioethics-Oct06.txt",
    "pg04"        : "pharmgkb-2004.txt",
    "pg07"        : "pharmgkb-070205.txt",
    "radiology"   : "daniel-radiology.txt",
    "mscannerbib" : "mscanner-bibliography.txt",
    "gdsmall"     : "genedrug-small.txt",
}

def query(*datasets):
    """Perform queries with per-feature pseudocount"""
    env = QueryEnvironment()
    rc.limit = 1000 # don't need a lot of results
    rc.threshold = 0
    rc.pseudocount = None # per-feature pseudocounts
    for dataset in datasets:
        rc.dataset = dataset
        env.standardQuery(rc.corpora / ds_map[dataset])
        
def retrieval(*datasets):
    """Carry out retrieval tests (adds -retrieval to the data set)"""
    env = QueryEnvironment()
    rc.limit = 1000
    rc.pseudocount = None # per-feature pseudocounts
    rc.retrieval_test_prop = 0.2
    for dataset in datasets:
        rc.dataset = dataset+"-retrieval"
        env.testRetrieval(rc.corpora / ds_map[dataset])

def heparin():
    """Performs a query using heparin (PubFinder example).
    Large inputs mask the effects of the choice of pseudocount,
    so I chose this small input to figure out what the best choice is:
    per-feature pseucocounts without a rarity cut-off.
    """
    env = QueryEnvironment()
    env.loadInput(rc.corpora / "pubfinder-heparin.txt")
    rc.limit = 500        # don't need a lot of results
    rc.threshold = 0      # lenient threshold of zero
    rc.pseudocount = None # per-feature pseudocounts
    for dataset, pseudocount, cutoff in [ 
        ("heparin-ps_const",    0.01, True),
        ("heparin-ps_per",      None, False),
        ("heparin-ps_constcut", 0.01, True),
        ("heparin-ps_percut",   None, True) ]:
        rc.dataset = dataset
        rc.pseudocount = pseudocount
        rc.cutoff = cutoff
        env.standardQuery()

def pharmdemo():
    """Special query which exports to the PharmDemo database for PharmGKB"""
    rc.dataset = "pharmdemo"
    rc.threshold = 50.0
    rc.limit = 10000
    env = QueryEnvironment()
    input = rc.corpora / "pharmgkb-070205.txt"
    input = rc.corpora / "genedrug-small.txt"
    env.getGDFilterResults(input, export_db=True)
        
if __name__ == "__main__":
    initLogger()
    if len(sys.argv) != 2:
        print "Please give python expression"
    else:
        eval(sys.argv[1])
