#!/usr/bin/python
import os
import sys
import SOAPpy

search_code = 9
http_proxy = os.environ.get("http_proxy", None)
server = SOAPpy.SOAPProxy("http://www.pharmgkb.org/services/SearchService", http_proxy=http_proxy)

result = cPickle.load(file("citations.pickle","rb"))
#result = server.specialSearch(SOAPpy.intType(search_code))

if len(result) == 0:
    raise "No results found!"
else:
    import cPickle
    #cPickle.dump(result, file("citations.pickle", "wb"), protocol=2)
    print result
    print len(result)
