"""Convert MEDLINE XML text/files into Article objects

@author: Graham Poulter
                                   

ArticleParser -- Convert MEDLINE XML text/files into Article objects

"""

import cPickle
import gzip
import logging as log
import pyRXPU
from article import Article

class ArticleParser:
    """Parse XML into article objects, optionally filtering MeSH terms"""

    def __init__(self, synonyms=None, exclude=None):
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
            if isinstance(synonyms, basestring):
                self.synonyms = cPickle.load(file(synonyms, "rb"))
            elif isinstance(synonyms, dict):
                self.synonyms = synonyms
            else:
                raise ValueError("Synonyms is neither a filename nor dictionary")
        if exclude is not None:
            if isinstance(exclude, basestring):
                self.exclude = cPickle.load(file(exclude, "rb"))
            elif isinstance(exclude, set):
                self.exclude = exclude
            else:
                raise ValueError("Excludes is neither a filename nor set")

    def __call__(self, text, check_id=None):
        """Call the L{parse} method"""
        return self.parse(text, check_id)

    def parse(self, text, check_id=None):
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
        NAME, ATTRS, CHILDREN, SPARE = range(0,4)
        parser = pyRXPU.Parser(Validate=0, ProcessDTD=0, TrustSDD=0)
        def getChildren(root, names):
            return (n for n in root[CHILDREN] if n[NAME] in names)
        def getChild(root, name):
            for n in root[CHILDREN]:
                if n[NAME] == name:
                    return n
        def parseCitation(root):
            if root[NAME] != 'MedlineCitation':
                raise ValueError("Attempted to parse non-MedlineCitation object")
            result = Article()
            for node1 in root[CHILDREN]:
                if node1[NAME] == 'PMID':
                    result.pmid = int(node1[CHILDREN][0])
                if node1[NAME] == 'Article':
                    for node2 in node1[CHILDREN]:
                        if node2[NAME] == 'ArticleTitle':
                            result.title = node2[CHILDREN][0]
                        if node2[NAME] == 'Abstract':
                            result.abstract = getChild(node2, 'AbstractText')[CHILDREN][0]
                        if node2[NAME] == 'Journal':
                            for node3 in node2[CHILDREN]:
                                if node3[NAME] == 'ISSN':
                                    result.issn = node3[CHILDREN][0]
                                if node3[NAME] == 'ISOAbbreviation':
                                    result.journal = node3[CHILDREN][0]
                                if node3[NAME] == 'JournalIssue':
                                    PubDate = getChild(node3, 'PubDate')
                                    if PubDate is None: continue
                                    Year = getChild(PubDate, 'Year')
                                    if Year is None: continue
                                    result.year = int(Year[CHILDREN][0])
                        if node2[NAME] == 'AuthorList':
                            for Author in getChildren(node2, ['Author']):
                                lastname = getChild(Author, 'LastName')
                                if lastname is None: continue
                                initials = getChild(Author, 'Initials')
                                if initials is None: continue
                                result.authors.append((initials[CHILDREN][0],lastname[CHILDREN][0]))
                if node1[NAME] == 'MeshHeadingList':
                    for MeshHeading in getChildren(node1, ['MeshHeading']):
                        for Descriptor in getChildren(
                            MeshHeading, ['DescriptorName', 'Descriptor', 'QualifierName', 'SubHeading']):
                            term = Descriptor[CHILDREN][0]
                            if term in self.synonyms:
                                log.debug("Detected synonym: %s --> %s",term,self.synonyms[term])
                                term = self.synonyms[term]
                            if term not in self.exclude:
                                result.meshterms.add(term)
                            else:
                                if term not in self.excludes_done:
                                    #self.log.debug("Detected exclude: %s",term)
                                    self.excludes_done.add(term)
                if node1[NAME] == 'ChemicalList':
                    for Chemical in getChildren(node1, ['Chemical']):
                        name = getChild(Chemical, 'NameOfSubstance')[CHILDREN][0]
                        result.chemicals.add(name)
            return result
        root = parser.parse(text)
        if root[NAME] == 'MedlineCitation' and root[ATTRS]['Status'] == 'MEDLINE':
            result = parseCitation(root)
            if check_id is not None:
                if result.pmid != check_id:
                    raise ValueError("Article PMID %s does not match check_id %s" % (result.pmid,check_id))
            yield result
        elif root[NAME] == 'PubmedArticle':
            for node1 in root[CHILDREN]:
                if node1[NAME] == 'MedlineCitation' and node1[ATTRS]['Status'] == 'MEDLINE':
                    result = parseCitation(node1)
                if check_id is not None:
                    if node1[NAME] == 'PubmedData':
                        matched = False
                        ArticleIdList = getChild(node1, 'ArticleIdList')
                        for ArticleId in getChildren(ArticleIdList, ['ArticleId']):
                            if ArticleId[CHILDREN][0] == check_id:
                                result.pmid = check_id
                                matched = True
                        if matched == False:
                            raise ValueError("Article IDs do not match check_id %s" % (check_id,))
            yield result
        elif root[NAME] == 'MedlineCitationSet':
            for MedlineCitation in getChildren(root, ['MedlineCitation']):
                if MedlineCitation[ATTRS]['Status'] == 'MEDLINE':
                    yield parseCitation(MedlineCitation)
        elif root[NAME] == 'PubmedArticleSet':
            for PubmedArticle in getChildren(root, ['PubmedArticle']):
                for MedlineCitation in getChildren(PubmedArticle, ['MedlineCitation']):
                    if MedlineCitation[ATTRS]['Status'] == 'MEDLINE':
                        yield parseCitation(MedlineCitation)
        else:
            raise ValueError("Could not find a valid PubMed citation")

    def parseFile(self, filename):
        """Return articles from an XML files

        @type filename: C{str}
        @param filename: .xml or .xml.gz MEDLINE citation archive

        @rtype: C{(Article)}
        @return: Parsed L{Article}'s
        """
        log.debug("Parsing XML file %s", filename.name)
        if filename.endswith(".gz"):
            text = gzip.open(filename, 'r').read()
        else:
            text = file(filename, 'r').read()
        for article in self.parse(text):
            #log.debug("Parsed article %d", article.pmid)
            yield article
