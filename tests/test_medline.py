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

from mscanner.medline.Article import Article
from mscanner.medline.FeatureDatabase import FeatureDatabase
from mscanner.medline.FeatureStream import FeatureStream
from mscanner.medline.FeatureMapping import FeatureMapping
from mscanner.medline.Updater import Updater
from mscanner import tests


class FeatureDatabaseTests(unittest.TestCase):
    
    def test_FeatureDatabase(self):
        d = FeatureDatabase()
        d.setitem(1, nx.array([1,3], d.ftype)) #eliminate duplicate features
        d.setitem(2, nx.array([2,3], d.ftype))
        self.assert_(all(d.getitem(1) == [1,3]))
        self.assert_(all(d.getitem(2) == [2,3]))
        self.assertRaises(KeyError, d.getitem, 3)
        self.failUnless(1 in d)
        self.failUnless(2 in d)
        self.failIf(3 in d)
        self.assertEqual(d.keys(), ['2','1'])
        self.assertEqual(list(d.__iter__()), ['2','1'])
        self.assertEqual(len(d), 2)
        d.delitem(2)
        self.failIf(2 in d)
        self.assertRaises(ValueError, d.setitem, 4, nx.array([3.3,4]))



class FeatureStreamTests(unittest.TestCase):

    def setUp(self):
        self.fn = path(tempfile.mktemp())
        
    def tearDown(self):
        if self.fn.isfile():
            self.fn.remove()
        
    def test_FeatureStream(self):
        with closing(FeatureStream(self.fn, nx.uint32, "a")) as fs:
            pmids = (12,34,56)
            dates = (20070101, 19980308, 20001207)
            feats = [nx.array([1,2,3,4], fs.ftype), 
                     nx.array([5,6,7,8], fs.ftype), 
                     nx.array([], fs.ftype)]
            for pmid, date, feat in zip(pmids, dates, feats):
                fs.additem(pmid, date, feat)
        with closing(FeatureStream(self.fn, nx.uint32, "r")) as fs:
            rpmids, rdates, rfeats = zip(*list(fs.iteritems()))
            self.assertEqual(pmids, rpmids)
            self.assertEqual(dates, rdates)
            for a, ra in zip(feats, rfeats):
                self.assert_(all(a == ra))



class FeatureMappingTests(unittest.TestCase):
    """For FeatureMapping class"""

    def setUp(self):
        self.fn = path(tempfile.mktemp())

    def tearDown(self):
        if self.fn.isfile():
            self.fn.remove()

    def test_FeatureMapping(self):
        fm = FeatureMapping(self.fn, nx.uint16)
        features = dict(Q=["A","B"], T=["A","C"])
        fm.add_article(features)
        self.assert_(nx.all(fm.get_vector(features) == [0,1,2,3]))
        self.assertEqual([fm[i] for i in [0,1,2,3,]], [("A","Q"), ("B","Q"),("A","T"),("C","T")])
        self.assertEqual(fm[1], ("B","Q"))
        self.assertEqual(fm[("C","T")], 3)
        self.assert_(nx.all(fm.counts == [1,1,1,1]))
        fm.dump()
        fm.load()
        fm.dump()
        fm.load()
        self.assertEqual(fm.features, [("A","Q"),("B","Q"),("A","T"),("C","T")])
        self.assertEqual(fm.feature_ids, {"Q":{"A":0,"B":1}, "T":{"A":2,"C":3}})
        self.assert_(nx.all(fm.class_mask("Q") == [1,1,0,0]))



class XMLParserTests(unittest.TestCase):

    def art_equal(self, a, b):
        for k, v in a.__dict__.iteritems():
            self.assertEqual(v, getattr(b, k))

    def test_parse_medline_xml(self):
        a1 = Article(
            pmid=1,
            title="T1",
            abstract="A1",
            journal="Mol. Biol. Rep.",
            issn="0301-4851",
            date_completed=(2000,6,29),
            pubyear=1999,
            meshterms=[("T1",),("T2",),("T3","Q4","Q5"),("T6","Q7")],
            authors=[("F1","L1"),("F2","L2")])
        a2 = Article(
            pmid=2,
            title="T2",
            abstract="A2",
            date_completed=(2000,6,29),
            meshterms=[("T1",),("T2",),("T3","Q4","Q5"),("T6","Q7")])
        b1, b2 = list(Article.parse_medline_xml(StringIO(xmltext)))
        self.art_equal(a1, b1)
        self.art_equal(a2, b2)



class UpdaterTests(unittest.TestCase):

    def setUp(self):
        self.home = path(tempfile.mkdtemp(prefix="medline-"))

    def tearDown(self):
        self.home.rmtree(ignore_errors=True)

    def test_Updater(self):
        from mscanner.configuration import rc
        h = self.home
        rc.cache = h
        xml = h/"test.xml"
        test_pmids = h/"pmids.txt"
        m = Updater.Defaults()
        xml.write_text(xmltext)
        m.add_directory(h, save_delay=1)
        logging.debug("".join((h/"articles.txt").lines()))
        test_pmids.write_lines(["1", "2"])
        a = m.adata.load_articles(test_pmids)
        logging.debug("Articles: %s", repr(a))
        logging.debug("%s", repr(m.fdata_all.featmap.features))
        self.assertEqual(a[0].pmid, 1)
        self.assertEqual(a[1].pmid, 2)
        self.assertEqual(m.fdata_mesh.featmap.counts, [2, 2, 2, 1, 2, 2, 2, 2])
        self.assertEqual(
            sorted(m.fdata_mesh.featmap.features), sorted([
                (u'Q4', 'qual'), (u'Q5', 'qual'), (u'Q7', 'qual'), 
                (u'T1', 'mesh'), (u'T2', 'mesh'), (u'T3', 'mesh'), (u'T6', 'mesh'), 
                (u'0301-4851', 'issn')]))



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
    unittest.main()
