#!/usr/bin/env python

"""Use SOAP to fetch current literature PMIDs from PharmGKB

Usage::
    python fetch-pharmgkb.py <outfile>

Prints a list of PubMed IDs retrieved from the PharmGKB web service,
and saves them to <outfile>.

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

import sys

def main(outfile):
    import os
    import SOAPpy
    search_code = 9
    http_proxy = os.environ.get("http_proxy", None)
    server = SOAPpy.SOAPProxy("http://www.pharmgkb.org/services/SearchService", http_proxy=http_proxy)
    result = server.specialSearch(SOAPpy.intType(search_code))
    if len(result) == 0:
        raise RuntimeError("No results found!")
    else:
        print "# %d PharmGKB PMIDs, going to %s" % (len(result), outfile)
        f = open(outfile,"w")
        for r in result:
            print r[7]
            f.write(r[7]+"\n")
        f.close()

if __name__ == "__main__":
    main(sys.argv[1])
