"""Test suite for mscanner.pharmdemo.Exporter

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from path import path
import sys
import tempfile
import unittest

from mscanner.medline.Article import Article
from mscanner.pharmdemo.Exporter import Exporter

articles=[
    Article(
    pmid=1000,
    title="TITLE1",
    abstract="ABSTRACT1",
    meshterms=["T1"],
    ),
    Article(
    pmid=2000,
    title="TITLE2",
    abstract="ABSTRACT2",
    meshterms=["T1"],
    )
]

articles[0].genedrug={'PKID1':["GENE1","GENE2"],'PKID2':["GENE2","GENE3"]}
articles[1].genedrug={'PKID1':["GENE2","GENE3"],'PKID2':["GENE3","GENE4"]}
        
class ExporterTests(unittest.TestCase):
    """Tests for Exporter module

    Tests: export_sqlfile
    Missing: exportSQlite, exportOracle*
    """
    def setUp(self):
        self.fn = path(tempfile.mktemp())

    def tearDown(self):
        self.fn.remove()
        
    def test(self):
        exporter = Exporter(articles);
        exporter.export_sqlfile(self.fn)
        """
        dbname = ":memory:"
        self.con = sqlite.connect(self.dbname)
        exporter.export_dbapi(self.con, articles)
        cur = self.con.cursor()
        cur.execute('SELECT * from genedrug order by id')

        rows=cur.fetchall()
        #for row in rows: print row
        assert rows==[
            (0, u'GENE3', u'PKID1', 1),
            (1, u'GENE2', u'PKID1', 2),
            (2, u'GENE2', u'PKID2', 1),
            (3, u'GENE4', u'PKID2', 1),
            (4, u'GENE1', u'PKID1', 1),
            (5, u'GENE3', u'PKID2', 2)]

        cur.execute('SELECT * from dg_pmids order by dg_id')
        rows = cur.fetchall()
        #for row in rows: print row
        assert rows==[
            (u'2000', None, 0),
            (u'1000', None, 1),
            (u'2000', None, 1),
            (u'1000', None, 2),
            (u'2000', None, 3),
            (u'1000', None, 4),
            (u'1000', None, 5),
            (u'2000', None, 5)]

        cur.execute('SELECT * from cbs order by pmid')
        rows=cur.fetchall()
        #for row in rows: print row
        assert rows==[
            (u'1000', u'TITLE1', u'ABSTRACT1', u'GENE1 GENE2 GENE3', u'PKID1 PKID2', None, None),
            (u'2000', u'TITLE2', u'ABSTRACT2', u'GENE2 GENE3 GENE4', u'PKID1 PKID2', None, None)]

        #self.con.close()
        if os.path.exists( self.dbname ):
            os.remove( self.dbname )
        """

if __name__ == "__main__":
    unittest.main()
