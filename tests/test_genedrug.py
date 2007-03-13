from path import path
import pprint 
import re
import tempfile
import unittest

from mscanner.article import Article
from mscanner import genedrug

pp = pprint.PrettyPrinter()

class ParseDrugsTest(unittest.TestCase):

    def test_parseDrugs(self):
        assert genedrug.parseDrugs(drugs_text) == parsed_drugs_correct

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

drugs_text = """PA10000\t17 beta-estradiol\t\t
PA10007\talbumin human\t\tAlbuminar-25|Albuminar-5|Albutein 25%|Albutein 5%|Buminate 25%|Buminate 5%|Plasbumin-25|Plasbumin-5|
PA10009\talefacept\t\tAmevive|
PA1001\t1-methyl-4-phenylpyridinium (MPP+)\t\t
PA10010\talemtuzumab\t\tCampath|
PA10011\talfacalcidol\t\tOne-Alpha|
PA10012\talteplase, recombinant\t\tActivase|Activase rt-PA|Cathflo Activase|"""

parsed_drugs_correct = {
    'PA10000': ['17 beta-estradiol'],
    'PA10007': ['albumin human',
                'Albuminar-25',
                'Albuminar-5',
                'Albutein 25%',
                'Albutein 5%',
                'Buminate 25%',
                'Buminate 5%',
                'Plasbumin-25',
                'Plasbumin-5'],
    'PA10009': ['alefacept', 'Amevive'],
    'PA1001': ['1-methyl-4-phenylpyridinium (MPP+)'],
    'PA10010': ['alemtuzumab', 'Campath'],
    'PA10011': ['alfacalcidol', 'One-Alpha'],
    'PA10012': ['alteplase, recombinant',
                'Activase',
                'Activase rt-PA',
                'Cathflo Activase']}

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
        self.genefinder = genedrug.CachingGeneFinder(self.fn)
        self.filter = genedrug.GeneDrugFilter({},self.genefinder)

    def tearDown(self):
        self.genefinder.__del__()
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

