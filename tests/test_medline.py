"""Test suite for mscanner.medline

                               

@license: This source file is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it and/or modify
it under the Do Whatever You Want Public License. Terms and conditions: 
   0. Do Whatever You Want
"""

from cStringIO import StringIO
from path import path
import tempfile
import unittest

from mscanner.article import Article
from mscanner.featuremap import FeatureMapping
from mscanner.medline import parse, MedlineCache, FileTracker
from mscanner.scorefile import getArticles
from mscanner.support.utils import usetempfile

class ArticleParserTests(unittest.TestCase):
    """
    Tests ArticleParser parse
    Not Testing: parseFile
    """
    def art_equal(self, a, b):
        self.assertEqual(a.pmid,b.pmid)
        self.assertEqual(a.title,b.title)
        self.assertEqual(a.abstract,b.abstract)
        self.assertEqual(a.journal,b.journal)
        self.assertEqual(a.issn,b.issn)
        self.assertEqual(a.year,b.year)
        self.assertEqual(a.meshterms,b.meshterms)
        self.assertEqual(a.authors,b.authors)
    def test(self):
        a1 = Article(1,"T1","A1","Mol. Biol. Rep.","0301-4851",1999,
                     [("T1",),("T2",),("T3","Q4","Q5"),("T6","Q7")],[("F1","L1"),("F2","L2")])
        a2 = Article(2,"T2","A2",None,None,None,
                     [("T1",),("T2",),("T3","Q4","Q5"),("T6","Q7")],[])
        result = list(parse(StringIO(xmltext)))
        self.art_equal(result[0],a1)
        self.art_equal(result[1],a2)


class MedlineCacheTests(unittest.TestCase):
    """
    Tests MedlineCache updateCacheFromDir (so also makeDBEnv, putArticleList)
    """
    def setUp( self ):
        self.home = path(tempfile.mkdtemp(prefix="medline-"))

    def tearDown( self ):
        self.home.rmtree(ignore_errors=True)

    def test( self ):
        h = self.home
        xml = h/"test.xml"
        pmids = h/"pmids.txt"
        artdb = h/"articles.db"
        featdb = h/"features.db"
        featstream = h/"features.stream"
        fmap = FeatureMapping(h/"featuremap.txt")
        m = MedlineCache(fmap,
                         h,
                         artdb,
                         featdb,
                         featstream,
                         h/"articles.txt",
                         h/"processed.txt",
                         h/"narticles.txt",
                         use_transactions=True,)
        xml.write_text(xmltext)
        m.updateCacheFromDir(h, save_delay=1)
        pmids.write_lines(["1", "2"])
        a = getArticles(artdb, pmids)
        print repr(a)
        self.assertEqual(a[0].pmid, 1)
        self.assertEqual(a[1].pmid, 2)
        self.assertEqual(fmap.counts, [2, 2, 2, 2, 2, 2, 2, 1])
        self.assertEqual(
            fmap.features, [
                (u'Q4', 'qual'), (u'Q5', 'qual'), (u'Q7', 'qual'), 
                (u'T1', 'mesh'), (u'T2', 'mesh'), (u'T3', 'mesh'), (u'T6', 'mesh'), 
                (u'0301-4851', 'issn')])
        
class FileTrackerTest(unittest.TestCase):

    @usetempfile
    def testFileTracker(self, fn):
        """For FileTracker.(__init__, add, toprocess, dump)"""
        t = FileTracker(fn)
        t.add(path("hack/a.xml"))
        t.add(path("cough/b.xml"))
        self.assertEqual(t.toprocess([path("foo/a.xml"), path("blah/c.xml")]), ["blah/c.xml"])
        t.dump()
        del t
        t = FileTracker(fn)
        self.assertEqual(t, set(['a.xml', 'b.xml']))
            
xmltext = u'''<?xml version="1.0"?>
<!DOCTYPE MedlineCitationSet PUBLIC "-//NLM//DTD Medline Citation, 1st January 2007//EN"
                                    "http://www.nlm.nih.gov/databases/dtd/nlmmedline_070101.dtd">
<MedlineCitationSet>

<MedlineCitation Owner="NLM" Status="MEDLINE">
<PMID>1</PMID>
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
<Article PubModel="Print">
<ArticleTitle>T2</ArticleTitle>
<Abstract>
<AbstractText>A2</AbstractText>
</Abstract>
</Article>
<ChemicalList>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance>C1</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>1</RegistryNumber>
<NameOfSubstance>C2</NameOfSubstance>
</Chemical>
</ChemicalList>
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
    unittest.main()
