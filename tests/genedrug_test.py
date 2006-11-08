"""Test genedrug module

@author: Graham Poulter
@contact: graham@cbio.uct.ac.za
                                   
"""

from path import path
import cPickle
import sys
import os
import unittest
import re
from pprint import PrettyPrinter

#http_proxy = "cache.uct.ac.za:8080"
http_proxy = None
pp=PrettyPrinter()
execfile("genedrug_data.py")

class GeneDrugTest(unittest.TestCase):

    def setUp(self):
        self.filter = genedrug.GeneDrugFilter(
            fulldrugs_correct,
            genedrug.CachingGeneFinder( genefinder_cache, http_proxy )
            )

    def testSentences(self):
        sentences=self.filter.listSentences(sentence_text)
        #print "SENTENCES: " + pp.pformat(sentences)
        assert sentences == sentences_correct

    def testDrugLoading(self):
        #print "FULLDRUGS: " + pp.pformat(self.filter.fulldrugs)
        #print "DRUGS: " + pp.pformat(self.filter.drugs)
        assert self.filter.fulldrugs == fulldrugs_correct
        assert self.filter.drugs == shortdrugs_correct
        
    def testGeneListing(self):
        genes=self.filter.listGenes( genedrug_text )
        #print "GENES: " + pp.pformat(genes)
        assert genes == genes_correct
        
    def testDrugListing(self):
        drugs = self.filter.listDrugs(drugs_text)
        #print "DRUGS: " + pp.pformat(drugs)
        assert drugs == drugs_correct

    def testGeneDrugListing(self):
        gd = self.filter.listGeneDrugs( genedrug_text )
        #print "GENEDRUG: " + pp.pformat(gd)
        #print "GENEFINDER_CACHE: " + pp.pformat( self.filter.geneFinder.cache )
        assert gd == genedrug_correct

    def testDrugListingThorough(self):
        """See what drugs we pull out from real articles."""
        # Get a caching gene-drug filter
        gdfilter=genedrug.GeneDrugFilter(
            load( file( "../var/drugs.cpickle", "rb" ) ),
            genedrug.CachingGeneFinder( "../var/gapscore.shelf", http_proxy)
            )
        """
        # Check that the sentences make sense
        for a in articles:
            print "Article %s" % (a.pmid,)
            sentences=gdfilter.listSentences(a.title+". "+a.abstract)
            sentencedrugs=gdfilter.listDrugs(s[0] for s in sentences)
            sen=0
            for mapping in sentencedrugs:
                #print "Sentence %d: %s\n"%(sen,sentences[sen][0])
                for (PKID,name) in mapping.iteritems():
                   print "Sentence %d: %s (%s) (%s)" % (sen,PKID,name,drugs[PKID][0])
                sen+=1
        """
        for a in articles:
            print "Article %s" % (a.pmid,)
            gdresult = gdfilter.listGeneDrugsArticle(a)
            print pp.pformat(gdresult)

def test_suite():
    suite=unittest.TestSuite()
    #suite.addTest(unittest.makeSuite(GeneDrugTest))
    suite.addTest(GeneDrugTest('testSentences'))
    suite.addTest(GeneDrugTest('testDrugLoading'))
    suite.addTest(GeneDrugTest('testDrugListing'))
    suite.addTest(GeneDrugTest('testGeneListing'))
    suite.addTest(GeneDrugTest('testGeneDrugListing'))
    #suite.addTest(GeneDrugTest('testDrugListingThorough'))
    return suite

if __name__=="__main__":
    unittest.main(defaultTest='test_suite')
