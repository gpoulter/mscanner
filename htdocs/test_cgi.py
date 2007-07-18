#!/export/home/medscan/local32/bin/python2.5
import web
import sys

sys.path.insert(0,"/export/home/medscan/source")
from mscanner import configuration

import web

urls = (
        '/(.*)', 'hello'
        )

class hello:
    def GET(self, name):
        i = web.input(times=1)
        if not name:
            name = 'world'
        for c in xrange(int(i.times)):
            print 'Hello,', name+'!'

if __name__ == "__main__":
    try:
        web.run(urls, globals())
    except KeyboardInterrupt:
        pass
