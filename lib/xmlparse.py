#!/usr/bin/env python

"""Convert MEDLINE XML text/files into Article objects

@author: Graham Poulter
                                   

ArticleParser -- Convert MEDLINE XML text/files into Article objects

"""

import cPickle
import gzip
import logging as log
import pyRXPU
import unittest
from article import Article

class ArticleParser:
    """Parse XML into article objects, optionally filtering MeSH terms"""

    def __init__( self, synonyms=None, exclude=None ):
        """Initialise parser
        
        @param synonyms: Name of a pickle mapping each MeSH term to a
        set of its synonyms.

        @param exclude: Name of a pickle containing a set of MeSH
        terms to exclude from articles while parsing.
        """
        self.synonyms = {}
        self.exclude = set()
        self.excludes_done = set()
        if synonyms is not None:
            if isinstance( synonyms, basestring ):
                self.synonyms = cPickle.load( file( synonyms, "rb" ) )
            elif isinstance( synonyms, dict ):
                self.synonyms = synonyms
            else:
                raise ValueError("Synonyms is neither a filename nor dictionary")
        if exclude is not None:
            if isinstance( exclude, basestring ):
                self.exclude = cPickle.load( file( exclude, "rb" ) )
            elif isinstance( exclude, set ):
                self.exclude = exclude
            else:
                raise ValueError("Excludes is neither a filename nor set")

    def __call__(self,text,check_id=None):
        """Call the L{parse} method"""
        return self.parse(text,check_id)

    def parse(self,text,check_id=None):
        """Yield articles from MEDLINE XML

        @type text: C{str}
        @param text: Text of one or more MedLine XML citations.
        
        @type check_id: C{str}
        @param check_id: A PubMed ID to check against article result.
        
        @raise ValueError: When citation could not be parsed, or does not
        have an identifier matching that supplied in check_id.
        
        @rtype: C{[Article]}
        @return: Article with pmid, title, abstract and meshterms fields
        set.
        """
        parser = pyRXPU.Parser(Validate=0,ProcessDTD=0,TrustSDD=0)
        NAME,ATTRS,CHILDREN,SPARE = range(0,4)
        def parseCitation(root):
            result = Article(pmid=0,title="",abstract="",meshterms=set())
            for node1 in root[CHILDREN]:
                if node1[NAME] == u'PMID':
                    result.pmid = int(node1[CHILDREN][0])
                if node1[NAME] == u'Article':
                    for node2 in node1[CHILDREN]:
                        if node2[NAME] == u'ArticleTitle':
                            result.title=node2[CHILDREN][0]
                        if node2[NAME] == u'Abstract':
                            for node3 in node2[CHILDREN]:
                                if node3[NAME] == u'AbstractText':
                                    result.abstract = node3[CHILDREN][0]
                if node1[NAME]==u'MeshHeadingList':
                    for MeshHeading in [n for n in node1[CHILDREN] if n[NAME] == 'MeshHeading']:
                        for termnode in [n for n in MeshHeading[CHILDREN] if n[NAME] in \
                        [u'DescriptorName',u'Descriptor',u'QualifierName',u'SubHeading'] ]:
                            term = termnode[CHILDREN][0]
                            if term in self.synonyms:
                                log.debug("Detected synonym: %s --> %s",term,self.synonyms[term])
                                term = self.synonyms[term]
                            if term not in self.exclude:
                                result.meshterms.add(term)
                            else:
                                if term not in self.excludes_done:
                                    #self.log.debug("Detected exclude: %s",term)
                                    self.excludes_done.add(term)
            return result
        root=parser.parse(text)
        if root[NAME]==u'MedlineCitation':
            result=parseCitation(root)
            if check_id is not None:
                if result.pmid!=check_id:
                    raise ValueError("Article PMID %s does not match check_id %s" % (result.pmid,check_id))
            yield result
        elif root[NAME]==u'PubmedArticle':
            for node1 in root[CHILDREN]:
                if node1[NAME]==u'MedlineCitation':
                    result=parseCitation(node1)
                if check_id is not None:
                    if node1[NAME] == u'PubmedData':
                        matched = False
                        ArticleIdList=[n for n in node1[CHILDREN] if n[NAME]=='ArticleIdList' ][0]
                        for ArticleId in [n for n in ArticleIdList[CHILDREN] if n[NAME]=='ArticleId']:
                            if ArticleId[CHILDREN][0]==check_id:
                                result.pmid=check_id
                                matched=True
                        if matched==False:
                            raise ValueError("Article IDs do not match check_id %s" % (check_id,))
            yield result
        elif root[NAME] == u'MedlineCitationSet':
            for MedlineCitation in [n for n in root[CHILDREN] if n[NAME]=='MedlineCitation']:
                yield parseCitation(MedlineCitation)
        else:
            raise ValueError("Could not find a valid PubMed citation")

    def parseFile( self, filename ):
        """Return articles from an XML files

        @type filename: C{str}
        @param filename: .xml or .xml.gz MEDLINE citation archive

        @rtype: C{(Article)}
        @return: Parsed L{Article}'s
        """
        log.debug( "Parsing XML file %s", filename.name )
        if filename.endswith( ".gz" ):
            text = gzip.open( filename, 'r' ).read()
        else:
            text = file( filename, 'r' ).read()
        for article in self.parse( text ):
            yield article

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

class _ArticleParserTests(unittest.TestCase):
    def art_equal( self, a, b):
        self.assertEqual(a.pmid,b.pmid)
        self.assertEqual(a.title,b.title)
        self.assertEqual(a.abstract,b.abstract)
        self.assertEqual(a.meshterms,b.meshterms)
    def test( self ):
        a1 = Article(1,"T","A",set(["T1","T2","T3","T4","T5","T6","T7"]))
        a2 = Article(2,"T","A",set(["T1","T2","T3","T4","T5","T6","T7"]))
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

if __name__ == "__main__":
    unittest.main()
