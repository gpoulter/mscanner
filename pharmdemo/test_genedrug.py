"""Test suite for mscanner.pharmdemo.genedrug

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from path import path
import operator
import pprint as pp
import re
import tempfile
import unittest

from mscanner.article import Article
from pharmdemo import genedrug

def log(string, *format):
    """Logging function, can avoid printing output"""
    print (string+"\n") % format


class DrugFinderTests(unittest.TestCase):
    """Tests for L{DrugFinder}"""

    def test_parse_drugs(self):
        drugs = genedrug.DrugFinder._parse_drugs(df_parser_input)
        log("DrugFinder._parse_drugs:\n %s", pp.pformat(drugs))
        self.assertEqual(drugs, df_parser_correct)

    def test_make_patterns(self):
        patterns = genedrug.DrugFinder._make_patterns(gd_drugtable)
        strpats = dict([ (n, r.pattern) for n, r in patterns.iteritems() ])
        log("DrugFinder._make_patterns:\n %s", strpats)
        self.assertEqual(strpats, df_patterns_correct)


class GeneFinderTests(unittest.TestCase):
    """Tests for L{GeneFinder}"""

    def setUp(self):
        self.fn = path(tempfile.mktemp(prefix="gd-"))
        self.genefinder = genedrug.GeneFinder(self.fn, 0.85)

    def tearDown(self):
        self.genefinder.close()
        self.fn.remove()

    def test_gapscore(self):
        genes = self.genefinder.gapscore(gf_text)
        log("GeneFinder.gapscore:\n %s", pp.pformat(genes))
        self.assertEqual(genes, gf_cache)
        self.assertEqual(self.genefinder.cache[repr(gf_text)], gf_cache)

    def test_find_genes(self):
        genes = self.genefinder.find_genes(gf_text)
        log("GeneFinder.find_genes:\n %s", pp.pformat(genes))
        self.assertEqual(genes, gf_correct)



class GeneDrugFinderTests(unittest.TestCase):
    """Tests for L{GeneDrugFinder}"""

    def setUp(self):
        gfinder = genedrug.GeneFinder(gd_genefinder_cache, 0.85)
        dfinder = genedrug.DrugFinder("")
        dfinder.drugtable = gd_drugtable
        dfinder.patterns = dfinder._make_patterns(dfinder.drugtable)
        self.gdfinder = genedrug.GeneDrugFinder({}, gfinder, dfinder)
        
    def tearDown(self):
        self.gdfinder.close()
        
    def test_find_drugs(self):
        drugs = self.gdfinder.drugfinder.find_drugs(gd_text)
        log("DrugFinder.find_drugs:\n %s", pp.pformat(drugs))
        self.assertEqual(drugs, gd_drugs_correct)

    def test_break_sentences(self):
        sentences = self.gdfinder.break_sentences(gd_text)
        log("GeneDrugFinder.break_sentences:\n %s", pp.pformat(sentences))
        self.assertEqual(sentences, gd_sentences_correct)

    def test_find_genedrugs_text(self):
        gd = self.gdfinder.find_genedrugs_text(gd_text)
        log("GeneDrugFinder.find_genedrugs_text:\n %s", pp.pformat(gd))
        self.assertEqual(gd, gd_genedrug_correct)


###############################################################################
### Data for DrugFinderTests

# Input for DrugFinder._parse_drugs
df_parser_input = """PA10000\t17 beta-estradiol\t\t
PA10007\talbumin human\t\tAlbuminar-25|Albuminar-5|Albutein 25%|Albutein 5%|Buminate 25%|Buminate 5%|Plasbumin-25|Plasbumin-5|
PA10009\talefacept\t\tAmevive|
PA1001\t1-methyl-4-phenylpyridinium (MPP+)\t\t
PA10010\talemtuzumab\t\tCampath|
PA10011\talfacalcidol\t\tOne-Alpha|
PA10012\talteplase, recombinant\t\tActivase|Activase rt-PA|Cathflo Activase|"""

# Output for DrugFinder._parse_drugs (pasted)
df_parser_correct = {
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

# Output for DrugFinder._make_patterns (pasted)
df_patterns_correct =  {
'D4': '\\b(gened)\\b',
'D2': '\\b(def\\s+xyz|bleh\\s+xyz)\\b',
'D1': '\\b(abcd|qirt)\\b',
}

###############################################################################
### Data for GeneFinderTests

# Input for GeneFinder.find_genes and gapscore
gf_text = """We observed an increase in mitogen-activated protein kinase (MAPK)
activity when administering Amevive. We found Hepatocyte Growth
Factorase (HGF) activity decreased when administering Campath."""

# End value of GeneFinder.cache (pasted)
gf_cache = [
            ('HGF',142,145,1),
            ('Hepatocyte Growth\nFactorase',113,140,1),
            ('MAPK',61,65,1),
            ('mitogen-activated protein kinase',27,59,1)
            ]

# Output for GeneFinder.find_genes (pasted)
gf_correct = [
    ('HGF', 142, 145),
    ('Hepatocyte Growth\nFactorase', 113, 140),
    ('MAPK', 61, 65), ('mitogen-activated protein kinase', 27, 59),
    ]

###############################################################################
### Data for GeneDrugFinder tests

# Genes present in gd_text
gd_genes = ["geneA", "geneB", "geneC", "geneD", "geneE"]

# Fill value for DrugFinder.drugtable
gd_drugtable = {
    "D1": ["Abc", "abc2", "abcd", "qirt"],
    "D2": ["def", "def xyz", "bleh xyz"],
    "D3": ["ghi", "ghi 25%"],
    "D4": ["geneD"],
}

# Sentences in gd_text
_sentences = [
    "One geneA Abc2. ",
    "One geneA Abcd2. ",
    "Dr. Smith geneB. ",
    "Two {geneC} qirt. ",
    "A 2.4 geneD? ",
    "An a.b geneE #bleh#@#$xyz# exclamation! ",
    "One. lowercase. ",
    "Just... an ellipsis.",
]

# Output for DrugFinder.find_drugs (pasted)
gd_drugs_correct =  [
    ('geneD', 74, 79, 'D4'),
    ('bleh#@#$xyz', 95, 106, 'D2'),
    ('Abcd', 26, 30, 'D1'),
    ('qirt', 62, 66, 'D1')
]

# Input for GeneDrugFinder.find_genedrugs_text
gd_text = "".join(_sentences)

# Output for GeneFinder.find_genes
gd_genes_correct = reduce(operator.__add__,
    [[(m.group(), m.start(), m.end()) for m in re.finditer(gene, gd_text)] \
    for gene in gd_genes])

# Output for GeneDrugFinder.find_genedrugs_text (pasted)
gd_genedrug_correct = {
    'def': set(['geneE']), 
    'Abc': set(['geneA', 'geneC']),
}

# Output for GeneDrugFinder.break_sentences
gd_sentences_correct = [ (s, gd_text.find(s), gd_text.find(s)+len(s)) for s in _sentences ]

# Fill value for GeneFinder.cache
gd_genefinder_cache = { repr(gd_text): [ (g,s,e,1.0) for g,s,e in gd_genes_correct ] }

###############################################################################

if __name__ == "__main__":
    unittest.main()

