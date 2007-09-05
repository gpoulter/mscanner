#!/usr/bin/env python

"""Query the MEDLINE database with a set of positive articles

Choose operation by executing a Python expression::

    python query.py 'query("pg04","pg07")'
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
import sys

from mscanner.configuration import rc, initLogger
from mscanner.queryenv import QueryEnvironment
            
ds_map = {
    "aids"        : "aids-bioethics-Oct06.txt",
    "pg04"        : "pharmgkb-2004.txt",
    "pg07"        : "pharmgkb-070205.txt",
    "radiology"   : "daniel-radiology.txt",
    "gdsmall"     : "genedrug-small.txt",
}
"""Mapping from dataset code to input file"""


def query(*datasets):
    """Perform queries with per-feature pseudocount"""
    env = QueryEnvironment()
    for dataset in datasets:
        rc.dataset = dataset
        env.standardQuery(rc.corpora / ds_map[dataset])


def retrieval(*datasets):
    """Carry out retrieval tests (adds -retrieval to the data set)"""
    env = QueryEnvironment()
    rc.limit = 1000 # Force 1000 results
    rc.retrieval_test_prop = 0.2
    for dataset in datasets:
        rc.dataset = dataset+"-retrieval"
        env.testRetrieval(rc.corpora / ds_map[dataset])


def pharmdemo():
    """Special query which exports to the PharmDemo database for PharmGKB"""
    rc.dataset = "pharmdemo"
    rc.threshold = 50.0 # Higher than usual threshold
    rc.limit = 10000 # Want lots of results
    env = QueryEnvironment()
    #input = rc.corpora / "pharmgkb-070205.txt"
    input = rc.corpora / "genedrug-small.txt"
    env.getGDFilterResults(input, export_db=True)


if __name__ == "__main__":
    initLogger()
    if len(sys.argv) != 2:
        print "Please provide a Python expression"
    else:
        eval(sys.argv[1])
