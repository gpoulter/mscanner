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
        ft = self.article.feats_iedb_concat()
        self.assertEqual(ft["iedb"], 
        ['DDD1', 'DDD2', 'DDD3', 'DDD6', 'AAA1.', 'LLL1', 
         'LLL2', 'QQQ1', 'TTT1', 'QQQ2', 'JJJ1', 'AAA2-FFF2'])
        logging.debug(str(ft))
        # Only use title/abstract for features
        ft = self.article.feats_iedb_word()
        self.assertEqual(ft["iedb"], 
        ['DDD1', 'AAA1.', 'TTT1', 'QQQ2', 'JJJ1', 'AAA2-FFF2'])

    def test_feats_mesh_qual_issn(self):
        ft = self.article.feats_mesh_qual_issn()
        logging.debug(str(ft))
        self.assertEqual(ft, {
            'qual': ['QQQ1', 'QQQ2'], 
            'issn': ['III1'], 
            'mesh': ['DDD1', 'DDD2', 'DDD3', 'DDD6']}
        )

    def test_feats_author(self):
        ft = self.article.feats_author()
        self.assertEqual(ft["a"], ['fff1 lll1', 'fff2 lll2'])
        logging.debug(str(ft))

    def test_feats_wmqia_filt(self):
        # Without filtering
        ft = self.article.feats_wmqia()
        self.assertEqual(ft["w"], 
          ['aaa1', 'ddd1', 'ttt1', 'qqq2', 'jjj1', 'aaa2-fff2'])
        logging.debug(str(ft))
        # With filtering
        ft = self.article.feats_wmqia_filt()
        self.assertEqual(ft["w"], 
          ['aaa1', 'ttt1', 'jjj1', 'aaa2-fff2'])
        logging.debug(str(ft))
        
    def test_word_extraction(self):
        """The other tests are just for consistencey. This one is to stress the
        feats_word tokenizer."""
        self.article.title = ""
        self.article.abstract = u"""EUROSENTINEL-study [ZCURVE adaptability].
        A/V-ATPase 85% +93.2% Nesb\xf8 ORF's 2,3-zn(ncs)2(c6h5nh2)2 stain--focus
        58.32 -a-3-5.3-b qe- te (RS) gene. M.G.T. D. radiodurans (orthologous)"""
        ft = self.article.feats_word()
        self.assertEqual(ft["w"], [u'gene', u'2,3-zn(ncs)2(c6h5nh2)2', u'rs',
        u'orthologous', u'radiodurans', u'a-3-5.3-b', u'a/v-atpase', u'zcurve',
        u'focus', u'eurosentinel-study', u'qe', u'nesb\xf8', u'orf',
        u'stain', u'te', u'adaptability', u'm.g.t'])
        logging.debug(str(ft))


if __name__ == "__main__":
    tests.start_logger()
    unittest.main()
