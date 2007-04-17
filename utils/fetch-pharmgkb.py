#!env python

"""Fetch current literature annotations from PharmGKB

Usage: python fetch-pharmgkb.py

Prints a list of PubMed IDs retrieved from the PharmGKB web service.
"""

import cPickle
import os
import sys
import SOAPpy

if __name__ == "__main__":
    
    search_code = 9
    http_proxy = os.environ.get("http_proxy", None)
    server = SOAPpy.SOAPProxy("http://www.pharmgkb.org/services/SearchService", http_proxy=http_proxy)
    result = server.specialSearch(SOAPpy.intType(search_code))
    
    if len(result) == 0:
        raise RuntimeError("No results found!")
    else:
        fname = "pharmgkb-pmids.txt"
        print "# %d PharmGKB PMIDs, going to %s" % (len(result), fname)
        f = file(fname,"w")
        for r in result:
            print r[7]
            f.write(r[7]+"\n")
        f.close()
            
