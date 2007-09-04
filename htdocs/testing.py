#!/export/home/medscan/local32/bin/python2.5
import web
import sys

sys.path.insert(0,"/export/home/medscan/source")
from mscanner import configuration

import web

urls = (
        '/hello/(.*)', 'hello',
        '/form', 'form',
        )

class hello:
    def GET(self, name):
        i = web.input(times=1)
        if not name:
            name = 'world'
        for c in xrange(int(i.times)):
            print 'Hello,', name+'!'

template = """
<html>
<head><title>Test Form</title</head>
<body>
<form action="/form" method="post">
<p>
<input type="checkbox" name="hidestuff">
<input type="submit">
</p>
</form>
</body>
</html>
"""
            
class form:
    def GET(self):
        print template
    def POST(self):
        import pprint as p
        p.pprint(web.input())

if __name__ == "__main__":
    try:
        web.run(urls, globals())
    except KeyboardInterrupt:
        pass
