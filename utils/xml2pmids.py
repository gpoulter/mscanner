import sys
sys.path.insert(0,"../../source/lib")
import xmlparse
from path import path
from cPickle import dump

filename = path( sys.argv[1] )
take = int(sys.argv[2])

outtxt = file( filename + ".txt", "w" )
parser = xmlparse.ArticleParser( "meshsynonyms.pickle", "meshexcludes.pickle"  )
articles = []
for i,x in enumerate( parser.parseFile( filename ) ):
    outtxt.write( "%d\n" % x.pmid )
    articles.append(x)
    if i >= take:
        break
dump( articles, file( filename + ".txt.pickle", "wb" ), protocol=2 )
