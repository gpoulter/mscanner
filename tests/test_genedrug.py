#!env python

from article import Article
import genedrug
import re
from path import path
import pprint 
import tempfile
import unittest

pp = pprint.PrettyPrinter()

class GeneDrugFilterTests(unittest.TestCase):
    """Tests for GeneDrugFilter class

    Tests: listSentences, listGenes, listDrugs, listGeneDrugs
    """

    def setUp(self):
        self.filter = genedrug.GeneDrugFilter(drugs,
        genedrug.CachingGeneFinder(genefinder_cache))

    def test_stripDrugs(self):
        drugtable = self.filter.stripDrugs(drugs)
        if(False):
            print "DRUGSTABLE: "
            for PKID, regex in drugtable.iteritems():
                print "%s: %s" % (PKID, regex.pattern)

    def test_listSentences(self):
        sentences = self.filter.listSentences(gdtext)
        #print "SENTENCES: " + pp.pformat(sentences)
        self.assertEqual(sentences, sentences_correct)

    def test_listGenes(self):
        genes = self.filter.listGenes(gdtext)
        #print "GENES: " + pp.pformat(genes)
        self.assertEqual(genes, genes_correct)

    def test_listDrugs(self):
        drugs = self.filter.listDrugs(gdtext)
        #print "DRUGS: " + pp.pformat(drugs)
        #self.assertEqual(drugs, drugs_correct)

    def test_listGeneDrugs(self):
        gd = self.filter.listGeneDrugs(gdtext)
        #print "GENEDRUG: " + pp.pformat(gd)
        self.assertEqual(gd, genedrug_correct)

drugs = {
    "D1": ["Abc", "abc2", "abcd", "qirt"],
    "D2": ["def", "def xyz", "bleh xyz"],
    "D3": ["ghi", "ghi 25%"],
    "D4": ["geneD"],
    }

genes = [
    "geneA",
    "geneB",
    "geneC",
    "geneD",
    "geneE",
    ]

sentences = [
    "One geneA Abc2. ",
    "One geneA Abcd2. ",
    "Dr. Smith geneB. ",
    "Two {geneC} qirt. ",
    "A 2.4 geneD? ",
    "An a.b geneE #bleh#@#$xyz# exclamation! ",
    "One. lowercase. ",
    "Just... an ellipsis.",
    ]

genedrug_correct = {'def': set(['geneE']), 'Abc': set(['geneA', 'geneC'])}

#### Auto-generated items

gdtext = "".join(sentences)

sentences_correct = [ (s, gdtext.find(s), gdtext.find(s)+len(s)) for s in sentences ]

genes_correct = []
for g in genes:
    genes_correct.extend( (m.group(), m.start(), m.end()) for m in re.finditer(g,gdtext) )

genefinder_cache = { repr(gdtext): [ (g,s,e,1.0) for g,s,e in genes_correct ] }

#####################################################################################

class GapScoreTests(unittest.TestCase):
    """Test the CachingGeneFinder using GapScore

    Tests: GeneDrugFilter.listGenes (implicitly testing CachingGeneFinder.findGenes)
    """

    def setUp(self):
        self.fn = path(tempfile.mktemp(prefix="gd-"))
        self.filter = genedrug.GeneDrugFilter({},genedrug.CachingGeneFinder(self.fn))

    def tearDown(self):
        del self.filter
        self.fn.remove()

    def test_listGenes(self):
        genes = self.filter.listGenes(gap_gdtext)
        #print "GENES: " + pp.pformat(genes)
        self.assertEqual(genes, gap_genes_correct)
        self.assertEqual(self.filter.geneFinder.cache[repr(gap_gdtext)], gap_genefinder_cache)
    
gap_gdtext="""We observed an increase in mitogen-activated protein kinase (MAPK)
activity when administering Amevive. We found Hepatocyte Growth
Factorase (HGF) activity decreased when administering Campath."""

gap_genefinder_cache = [
            ('HGF',142,145,1),
            ('Hepatocyte Growth\nFactorase',113,140,1),
            ('MAPK',61,65,1),
            ('mitogen-activated protein kinase',27,59,1)
            ]

gap_genes_correct = [
    ('HGF', 142, 145),
    ('Hepatocyte Growth\nFactorase', 113, 140),
    ('MAPK', 61, 65), ('mitogen-activated protein kinase', 27, 59),
    ]

if __name__ == "__main__":
    unittest.main()
