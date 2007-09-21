#!/usr/bin/env python

"""Performs queries using datasets from the MScanner paper

Example::

    python query.py pg04 pg07
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

import sys

from mscanner.configuration import rc, initLogger
from mscanner import scorefile, queryenv


dataset_map = {
    "aids"        : "aids-bioethics-Oct06.txt",
    "pg04"        : "pharmgkb-2004.txt",
    "pg07"        : "pharmgkb-070205.txt",
    "radiology"   : "daniel-radiology.txt",
    "gdsmall"     : "genedrug-small.txt",
    "invalid"     : "testing-invalid.txt",
}
"""Mapping from dataset code to input file"""


def do_query(*datasets):
    env = scorefile.Databases()
    for dataset in datasets:
        if dataset not in dataset_map:
            raise ValueError("Invalid Data Set %s" % dataset)
        rc.dataset = dataset
        query = queryenv.Query(rc.working / "query" / rc.dataset, env)
        query.query(rc.corpora / dataset_map[dataset])
    env.close()


if __name__ == "__main__":
    initLogger()
    rc.citations_per_file = 100
    do_query(*sys.argv[1:])
