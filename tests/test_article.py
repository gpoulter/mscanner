"""Test suite for mscanner.medline.Article

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from __future__ import with_statement
from contextlib import closing
import logging
import unittest

from mscanner import tests
from mscanner.medline.Article import Article


class ArticleTests(unittest.TestCase):
    """Test feature extraction from articles"""
    
    def setUp(self):
        self.article = Article(
            pmid=1,
            title="TTT1 DDD1",
            abstract="AAA1. AAA2-FFF2 QQQ2 JJJ1 TTT1",
            journal="JJJ1",
            issn="III1",
            date_completed=(2000,6,29),
            pubyear=1999,
            meshterms=[("DDD1",),("DDD2",),("DDD3","QQQ1","QQQ2"),("DDD6","QQQ2")],
            authors=[("FFF1","LLL1"),("FFF2","LLL2")])
    
    def test_feats_iedb(self):
        ft = sorted(self.article.feats_iedb_concat()["w"])
        self.assertEqual(ft, 
            ['AAA1.', 'AAA2-FFF2', 'DDD1', 'DDD2', 'DDD3',
            'DDD6', 'JJJ1', 'LLL1', 'LLL2', 'QQQ1', 'QQQ2', 'TTT1'])
        logging.debug(str(ft))
        # Only use title/abstract for features
        ft = sorted(self.article.feats_iedb_word()["w"])
        self.assertEqual(ft, 
            ['AAA1.', 'AAA2-FFF2', 'DDD1', 'JJJ1', 'QQQ2', 'TTT1'])
        
    def test_feats_mesh_qual_issn(self):
        ft = self.article.feats_mesh_qual_issn()
        logging.debug(str(ft))
        self.assertEqual(ft, {
            'qual': ['QQQ1', 'QQQ2'], 
            'issn': ['III1'], 
            'mesh': ['DDD1', 'DDD2', 'DDD3', 'DDD6']}
        )

    def test_feats_author(self):
        ft = sorted(self.article.feats_author()["a"])
        self.assertEqual(ft, ['fff1 lll1', 'fff2 lll2'])
        logging.debug(str(ft))

    def test_feats_wmqia_filt(self):
        # Without filtering
        ft = sorted(self.article.feats_wmqia()["w"])
        self.assertEqual(ft,
        ['AAA1', 'AAA2-FFF2', 'DDD1', 'JJJ1', 'QQQ2', 'TTT1'])
        logging.debug(str(ft))
        # With filtering
        ft = sorted(self.article.feats_wmqia_filt()["w"])
        self.assertEqual(ft, 
        ['AAA1', 'AAA2-FFF2', 'JJJ1', 'TTT1'])
        logging.debug(str(ft))
        
        
class WordExtractionTests(unittest.TestCase):
    """Test ways of extracting words from text"""

    def setUp(self):
        self.article = Article(pmid=1, title="", abstract= u"""335-349 EU-Ses
        [ZC ad]. A/V-ATP +85.32% Neb\xf8 n=>3 -ORF's 2,3-zn(ncs)2 st--foc :-Foc:
        -a-3.-b qe-- (RS) M.G.T. (orth)""")

    def test_word(self):
        self.assertEqual(sorted(self.article.feats_word()["w"]),
        [u'2,3-zn(ncs)2', u'A-3.-b', u'Ad', u'EU-Ses', u'Foc',
         u'M.G.T', u'Neb\xf8', u'ORF', u'Orth', u'Qe', u'RS', u'St', u'V-ATP', u'ZC'])
        
    def test_word_fold(self):
        self.assertEqual(sorted(self.article.feats_word_fold()["w"]),
        [u'2,3-zn(ncs)2', u'a-3.-b',  u'ad', u'eu-ses', u'foc',
         u'm.g.t', u'neb\xf8', u'orf', u'orth', u'qe', u'rs', u'st', u'v-atp', u'zc'])
        
    def test_word_num(self):
        self.assertEqual(sorted(self.article.feats_word_num()["w"]),
        [u'2,3-zn(ncs)2', u'335-349', u'85.32', u'A-3.-b', u'Ad', u'EU-Ses',
         u'Foc', u'M.G.T', u'Neb\xf8', u'ORF', u'Orth', u'Qe', u'RS', u'St',
         u'V-ATP', u'ZC'])

    def test_word_strip(self):
        self.assertEqual(sorted(self.article.feats_word_strip()["w"]), [u'32',
        u'335', u'349', u'85', u'ad', u'atp', u'eu', u'foc', u'ncs',
        u'neb\xf8', u'orf', u'orth', u'qe', u'rs', u'ses', u'st', u'zc',
        u'zn'])

    def test_iedb_word(self):
        self.assertEqual(sorted(self.article.feats_iedb_word()["w"]), [u'(RS)',
            u'(orth)', u'+85.32%', u"-ORF's", u'-a-3.-b', u'2,3-zn(ncs)2',
            u'335-349', u':-Foc:', u'A/V-ATP', u'EU-Ses', u'M.G.T.', u'Neb\xf8',
            u'[ZC', u'ad].',  u'n=>3', u'qe--', u'st--foc'])


if __name__ == "__main__":
    tests.start_logger()
    unittest.main()
