"""Test suite for mscanner.medline

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from __future__ import with_statement
from contextlib import closing

from cStringIO import StringIO
import logging
import numpy as nx
from path import path
import tempfile
import unittest

from mscanner.configuration import rc
from mscanner.medline.Article import Article
from mscanner.medline.FeatureData import FeatureData
from mscanner.medline.FeatureVectors import FeatureVectors
from mscanner.medline.FeatureStream import FeatureStream
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.medline.Updater import Updater
from mscanner.scripts import update
from mscanner import tests


class FeatureVectorsTests(unittest.TestCase):
    """FeatureVector SQLite database tests"""
    
    def test(self):
        a1 = (1, 19990101, [1,3])
        a2 = (2, 20000101, [2,3,5434,3243232])
        a3 = (3, 20010101, [5,10,15])
        a2b = (2, 20020101, [5,6,7])
        # Test adding 1,2,3
        d = FeatureVectors(filename=None)
        for a in a1,a2,a3: d.add_record(*a) 
        self.failUnless(d.get_records([1]).next() == a1)
        self.failUnless(d.get_records([2]).next() == a2)
        self.assertEqual(list(d.get_records([4])), [])
        self.failUnless(1 in d)
        self.failUnless(2 in d)
        self.failIf(4 in d)
        self.assertEqual(len(d), 3)
        # Test iteration
        filtered = list(d.iteritems(20000101,20020101))
        self.assertEqual(len(filtered), 2)
        self.assertEqual(filtered[0], a2)
        # Test updating
        d.update_record(*a2b)
        self.failUnless(d.get_records([2]).next() == a2b)
        # Test random-get
        self.assertEqual(sorted(d.get_random(3, 20000101)), [a2b,a3])
        # Test multi-get
        self.assertEqual(list(d.get_records([3,2])), [a2b,a3])
        

class FeatureStreamTests(unittest.TestCase):

    def setUp(self):
        self.fn = path(tempfile.mktemp())
        
    def tearDown(self):
        if self.fn.isfile():
            self.fn.remove()
        
    def test(self):
        """Write to FeatureStream, and read it all back."""
        pmids = (12,34,56)
        dates = (20070101, 19980308, 20001207)
        feats = [[1,2,6,5484], [5,6,8342,9000,9001], []]
        with closing(FeatureStream(self.fn, rdonly=False)) as fs:
            for pmid, date, feat in zip(pmids, dates, feats):
                fs.additem(pmid, date, feat)
        with closing(FeatureStream(self.fn, rdonly=True)) as fs:
            rpmids, rdates, rfeats = zip(*list(fs.iteritems()))
            self.assertEqual(pmids, rpmids)
            self.assertEqual(dates, rdates)
            for a, ra in zip(feats, rfeats):
                self.assertEqual(a, ra)



class FeatureMappingTests(unittest.TestCase):

    def test(self):
        """FeatureMapping - store and retrieve instances."""
        fm = FeatureMapping(filename=None)
        a1 = dict(Q=["A","B"], T=["A","C"])
        # Test adding of a feature vector
        fm.add_article(a1)
        self.assertEqual(fm.make_vector(a1), [1,2,3,4])
        self.assertEqual([fm.get_feature(i) for i in [1,2,3,4]], [("A","Q"), ("B","Q"),("A","T"),("C","T")])
        self.assert_(nx.all(fm.counts == [0,1,1,1,1]))
        self.assert_(nx.all(fm.type_mask("Q") == [1,1,1,0,0]))
        # Test read/write of unicode characters
        a2 = {"Q":["B", u"D\xd8"]}
        fm.add_article(a2)
        logging.debug("%s", str(fm.get_feature(5)))
        # Test removing of a feature vector
        self.assert_(nx.all(fm.counts == [0,1,2,1,1,1]))
        fm.remove_vector(fm.make_vector(a2))
        del fm._counts
        self.assert_(nx.all(fm.counts == [0,1,1,1,1,0]))
        fm.con.close()



class XMLParserTests(unittest.TestCase):

    def art_equal(self, a, b):
        for k, v in a.__dict__.iteritems():
            self.assertEqual(v, getattr(b, k))

    def test_parse_medline_xml(self):
        """Parse Medline XML and compare against correct Article object."""
        a_correct = Article(
            pmid=1,
            title="T1",
            abstract="A1",
            journal="Mol. Biol. Rep.",
            issn="0301-4851",
            date_completed=(2000,6,29),
            pubyear=1999,
            meshterms=[("T1",),("T2",),("T3","Q4","Q5"),("T6","Q7")],
            authors=[("F1","L1"),("F2","L2")])
        b_correct = Article(
            pmid=2,
            title="T2",
            abstract="A2",
            date_completed=(2000,6,29),
            meshterms=[("T1",),("T2",),("T3","Q4","Q5"),("T6","Q7")])
        a, b = list(Article.parse_medline_xml(StringIO(xmltext)))
        self.art_equal(a, a_correct)
        self.art_equal(b, b_correct)


class FeatureDataTests(unittest.TestCase):
    """Tests of FeatureData"""
    
    def setUp(self):
        self.home = path(tempfile.mkdtemp(prefix="fdata-"))
        rc.articles_home = self.home

    def tearDown(self):
        self.home.rmtree(ignore_errors=True)
        
    def test(self):
        """Tests of FeatureData"""
        articles = {
            "333": Article(333,date_completed=(1990,01,01),meshterms=[("A","B"),"C"]),
            "444": Article(444,date_completed=(1990,01,01),meshterms=[("D","B"),"E"]),
            }
        fd = FeatureData.Defaults("feats_mesh_qual_issn", False)
        fd.add_articles(articles.itervalues())
        fd.close()
        fd.featuredb.filename.remove()
        fd.featmap.filename.remove()
        fd.fstream.filename.remove()
        fd = FeatureData.Defaults("feats_mesh_qual_issn", False)
        fd.regenerate(articles)
        self.assertEqual(len(fd.featmap), 6)
        self.assertEqual(len(fd.featuredb), 2)
        fd.close()


class UpdaterTests(unittest.TestCase):
    """Tests L{Updater}"""

    def setUp(self):
        self.home = path(tempfile.mkdtemp(prefix="medline-"))
        rc.articles_home = self.home
        rc.medline = self.home

    def tearDown(self):
        self.home.rmtree(ignore_errors=True)

    def test(self):
        """Parse a sample directory and compare against known feature counts.
        """
        h = self.home
        update.featurespaces = ["feats_mesh_qual_issn", "feats_wmqia"]
        (h/"test.xml").write_text(xmltext)
        # Tests Updater.Defaults and add_directory
        m = update.medline()
        a = [ m.artdb[x] for x in ["1","2"] ]
        #logging.debug("Articles: %s", repr(a))
        self.assertEqual(a[0].pmid, 1)
        self.assertEqual(a[1].pmid, 2)
        self.assert_(nx.all(m.fdata_list[0].featmap.counts == [0, 2, 2, 2, 1, 2, 2, 2, 2]))
        self.assertEqual(
            sorted(m.fdata_list[0].featmap.con.execute("SELECT name,type FROM fmap").fetchall()), 
                sorted([(u'', u''),
                (u'Q4', 'qual'), (u'Q5', 'qual'), (u'Q7', 'qual'), 
                (u'T1', 'mesh'), (u'T2', 'mesh'), (u'T3', 'mesh'), (u'T6', 'mesh'), 
                (u'0301-4851', 'issn')]))
        m.close()
        # Test regeneration
        logging.debug("Removing MeSH featuredb")
        m.fdata_list[0].featuredb.filename.remove()
        m = update.regenerate()
        self.assertEqual(len(m.fdata_list[0].featuredb), 2)
        m.close()
        # Test save/load pickle
        (h/"pmids.txt").write_lines(["1","2"])
        update.save_pickle(h/"pmids.txt")
        m = update.load_pickles(h/"pmids.pickle.gz")



xmltext = u'''<?xml version="1.0"?>
<!DOCTYPE MedlineCitationSet PUBLIC 
"-//NLM//DTD Medline Citation, 1st January 2007//EN"
"http://www.nlm.nih.gov/databases/dtd/nlmmedline_070101.dtd">
<MedlineCitationSet>

<MedlineCitation Owner="NLM" Status="MEDLINE">
<PMID>1</PMID>
<DateCompleted>
 <Year>2000</Year><DateCompleted>
    <Year>2000</Year>
    <Month>06</Month>
    <Day>29</Day>
</DateCompleted>

 <Month>06</Month>
 <Day>29</Day>
</DateCompleted>
<Article PubModel="Print">
 <Journal>
 <ISSN IssnType="Print">0301-4851</ISSN>
 <JournalIssue CitedMedium="Print">
 <Volume>26</Volume>
 <Issue>3</Issue>
 <PubDate>
 <Year>1999</Year>
 <Month>Aug</Month>
 </PubDate>
 </JournalIssue>
 <Title>Molecular biology reports. </Title>
 <ISOAbbreviation>Mol. Biol. Rep.</ISOAbbreviation>
 </Journal>
 <ArticleTitle>T1</ArticleTitle>
 <Abstract>
 <AbstractText>A1</AbstractText>
 </Abstract>
 <AuthorList CompleteYN="Y">
 <Author ValidYN="Y">
  <LastName>L1</LastName>
  <ForeName>T J</ForeName>
 <Initials>F1</Initials>
 </Author>
 <Author ValidYN="Y">
  <LastName>L2</LastName>
  <ForeName>A J</ForeName>
 <Initials>F2</Initials>
 </Author>
 </AuthorList>
</Article>
<MedlineJournalInfo>
<Country>ENGLAND</Country>
<MedlineTA>Mol. Biol. Rep.</MedlineTA>
<NlmUniqueID>8712028</NlmUniqueID>
</MedlineJournalInfo>
<MeshHeadingList>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T1</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="Y">T2</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T3</DescriptorName>
<QualifierName MajorTopicYN="N">Q4</QualifierName>
<QualifierName MajorTopicYN="Y">Q5</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T6</DescriptorName>
<QualifierName MajorTopicYN="N">Q7</QualifierName>
</MeshHeading>
</MeshHeadingList>
</MedlineCitation>

<MedlineCitation Owner="NLM" Status="MEDLINE">
<PMID>2</PMID>
<DateCompleted>
 <Year>2000</Year>
 <Month>06</Month>
 <Day>29</Day>
</DateCompleted>
<Article PubModel="Print">
 <ArticleTitle>T2</ArticleTitle>
 <Abstract>
 <AbstractText>A2</AbstractText>
 </Abstract>
</Article>
<MeshHeadingList>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T1</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="Y">T2</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T3</DescriptorName>
<QualifierName MajorTopicYN="N">Q4</QualifierName>
<QualifierName MajorTopicYN="Y">Q5</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T6</DescriptorName>
<QualifierName MajorTopicYN="N">Q7</QualifierName>
</MeshHeading>
</MeshHeadingList>
</MedlineCitation>

</MedlineCitationSet>
'''

if __name__ == "__main__":
    tests.start_logger()
    #unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(ArticleTests)
    #suite = unittest.TestLoader().loadTestsFromTestCase(FeatureMappingTests)
    #suite = unittest.TestLoader().loadTestsFromTestCase(FeatureVectorsTests)
    suite = unittest.TestLoader().loadTestsFromTestCase(UpdaterTests)
    unittest.TextTestRunner(verbosity=2).run(suite)


