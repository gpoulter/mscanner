#!env python

from article import Article
from path import path
import unittest
from xmlparse import ArticleParser

class ArticleParserTests(unittest.TestCase):
    """
    Tests ArticleParser parse
    Not Testing: parseFile
    """
    def art_equal( self, a, b):
        self.assertEqual(a.pmid,b.pmid)
        self.assertEqual(a.title,b.title)
        self.assertEqual(a.abstract,b.abstract)
        self.assertEqual(a.meshterms,b.meshterms)
        self.assertEqual(a.chemicals,b.chemicals)
    def test( self ):
        a1 = Article(1,"T","A",set(["T1","T2","T3","T4","T5","T6","T7"]),set([]))
        a2 = Article(2,"T","A",set(["T1","T2","T3","T4","T5","T6","T7"]),set(["C1","C2"]))
        parser = ArticleParser()
        result = list( parser.parse( xmltext ) )
        self.art_equal(result[0],a1)
        self.art_equal(result[1],a2)
        synonyms = {"T1":"T2"}
        exclude = set(["T3","T4"])
        parser = ArticleParser( synonyms, exclude )
        result = list( parser.parse( xmltext ) )
        a1.meshterms.remove("T1")
        a1.meshterms.remove("T3")
        a1.meshterms.remove("T4")
        self.art_equal(result[0],a1)

xmltext = u'''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE MedlineCitationSet PUBLIC "-//NLM//DTD Medline Citation, 1st January 2006//EN"
                                    "http://www.nlm.nih.gov/databases/dtd/nlmmedline_060101.dtd">
<MedlineCitationSet>

<MedlineCitation Owner="NLM" Status="MEDLINE">
<PMID>1</PMID>
<Article PubModel="Print">
<ArticleTitle>T</ArticleTitle>
<Abstract>
<AbstractText>A</AbstractText>
</Abstract>
</Article>
<MeshHeadingList>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T1</DescriptorName>
</MeshHeading>
<MeshHeading>
<Descriptor MajorTopicYN="Y">T2</Descriptor>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T3</DescriptorName>
<QualifierName MajorTopicYN="N">T4</QualifierName>
<QualifierName MajorTopicYN="Y">T5</QualifierName>
</MeshHeading>
<MeshHeading>
<Descriptor MajorTopicYN="N">T6</Descriptor>
<SubHeading MajorTopicYN="N">T7</SubHeading>
</MeshHeading>
</MeshHeadingList>
</MedlineCitation>

<MedlineCitation Owner="NLM" Status="MEDLINE">
<PMID>2</PMID>
<Article PubModel="Print">
<ArticleTitle>T</ArticleTitle>
<Abstract>
<AbstractText>A</AbstractText>
</Abstract>
</Article>
<ChemicalList>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance>C1</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance>C2</NameOfSubstance>
</Chemical>
</ChemicalList>
<MeshHeadingList>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T1</DescriptorName>
</MeshHeading>
<MeshHeading>
<Descriptor MajorTopicYN="Y">T2</Descriptor>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">T3</DescriptorName>
<QualifierName MajorTopicYN="N">T4</QualifierName>
<QualifierName MajorTopicYN="Y">T5</QualifierName>
</MeshHeading>
<MeshHeading>
<Descriptor MajorTopicYN="N">T6</Descriptor>
<SubHeading MajorTopicYN="N">T7</SubHeading>
</MeshHeading>
</MeshHeadingList>
</MedlineCitation>

</MedlineCitationSet>
'''

if __name__ == "__main__":
    unittest.main()
