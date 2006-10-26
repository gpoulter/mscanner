#!env python

"""Fetch current literature annotations from PharmGKB

Usage: python fetch-pharmgkb-pmids.py

Prints a list of PubMed IDs retrieved from the PharmGKB web service.
"""

import cPickle
import os
import sys
import SOAPpy

search_code = 9
http_proxy = os.environ.get("http_proxy", None)
server = SOAPpy.SOAPProxy("http://www.pharmgkb.org/services/SearchService", http_proxy=http_proxy)
result = server.specialSearch(SOAPpy.intType(search_code))

if len(result) == 0:
    raise RuntimeError("No results found!")
else:
    print "# %s PharmGKB evidence PubMed IDs" % len(result)
    for r in result:
        print r[7]
