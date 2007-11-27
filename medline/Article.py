"""Provides the Article class"""


                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


class Article:
    """Database record for a Medline citation.
    
    The article is converted to a tuple which stored in a Berkeley DB
    indexed by PubMed ID.

    @ivar pmid: PubMed ID of the article (int)

    @ivar title: Title of the article (string)

    @ivar abstract: Abstract of the article (sring)

    @ivar journal: Medline abbreviated journal title (string)

    @ivar issn: Journal ISSN code (string)

    @ivar date_completed: (year,month,day) as integers
    
    @ivar year: Year of publication (int)
    
    @ivar meshterms: MeSH as a list of (descriptor, qual, qual, ...) tuples

    @ivar authors: Authors as a list of (initials, lastname) tuples of 
    """
    
    def __init__(self,
                 pmid=None,
                 title=None,
                 abstract=None,
                 journal=None,
                 issn=None,
                 date_completed=None,
                 pubyear=None,
                 meshterms=None,
                 authors=None):
        """Constructor, where parameters set instance variables."""
        self.pmid = pmid
        self.title = title
        self.abstract = abstract
        self.journal = journal
        self.issn = issn
        self.date_completed = date_completed
        self.pubyear = pubyear
        self.meshterms = meshterms if meshterms else []
        self.authors = authors if authors else []


    def __repr__(self):
        """Evaluatable representation of the object"""
        import pprint as pp
        s = ["Article("]
        for k,v in self.__dict__.iteritems():
            s.append("%s=%s" % (k, pp.pformat(v)))
        s[-1] += ")"
        return "\n".join(s)


    @staticmethod
    def parse_medline_xml(stream):
        """Generate Article objects by parsing a Medline XML file
        
        @param stream: File-like object of MedlineCitation XML
        
        @return: Iteratation over parsed Article objects
        """
        import xml.etree.cElementTree as ET
        context = ET.iterparse(stream, events=("start", "end"))
        context = iter(context)
        event, root = context.next()
        for event, record in context:
            if event == "end" and record.tag == "MedlineCitation":
                if record.get("Status") == "MEDLINE":
                    r = Article()
                    r.pmid = int(record.findtext("PMID"))
                    dc = record.find("DateCompleted")
                    r.date_completed = (
                        int(dc.findtext("Year")),
                        int(dc.findtext("Month")),
                        int(dc.findtext("Day")))
                    art = record.find("Article")
                    r.issn = art.findtext("Journal/ISSN")
                    r.pubyear = art.findtext("Journal/JournalIssue/PubDate/Year")
                    if r.pubyear is not None:
                        r.pubyear = int(r.pubyear)
                    r.title = art.findtext("ArticleTitle")
                    r.abstract = art.findtext("Abstract/AbstractText")
                    r.authors = [(a.findtext("Initials"), a.findtext("LastName")) 
                                      for a in art.findall("AuthorList/Author")]
                    r.journal = record.findtext("MedlineJournalInfo/MedlineTA")
                    for heading in record.findall("MeshHeadingList/MeshHeading"):
                        descriptor = heading.findtext("DescriptorName")
                        quals = [ q.text for q in heading.findall("QualifierName") ]
                        r.meshterms.append(tuple([descriptor] + quals))
                    yield r
                root.clear()


    def mesh_features(self):
        """Return MeSH and Journal features in this article, as a mapping from
        feature type to feature strings of that type. Types are 'mesh', 'qual'
        and 'issn'."""
        # Get MeSH headings, qualifiers and ISSN from article
        headings, quals = [], []
        for term in self.meshterms:
            headings.append(term[0])
            if len(term) > 1:
                for q in term[1:]:
                    if q not in quals:
                        quals.append(q)
        issns = [self.issn] if self.issn is not None else []
        # Get the feature vector while possibly them to the feature mapping
        return dict(mesh=headings, qual=quals, issn=issns)


    def word_features(self, stopwords):
        """Return a list of word features present in the article."""
        import re
        text = self.title.lower() + " "
        if self.abstract is not None:
            text += self.abstract.lower()
        # Get rid of non-alphabetics, and multiple spaces
        losechars = re.compile(r'[^a-z]+')
        text = losechars.sub(' ', text).strip()
        text = re.sub(r'\s+', ' ', text)
        # Keep only longish words not in stopwords
        return [x for x in text.split() if len(x) >= 3 and x not in stopwords]


    def all_features(self, stopwords):
        """Adds L{word_features} to L{mesh_features} using the 'word' feature
        type."""
        features = self.mesh_features()
        features["word"] = self.word_features(stopwords)
        return features

    