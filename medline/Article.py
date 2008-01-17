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

import logging
import re

# List of words to ignore in free text
stopwords = []

def init_stopwords(filename):
    """Initialise the L{stopwords} module variable."""
    global stopwords
    if filename is not None:
        stopwords = filename.lines(retain=False)


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


    def feats_mesh_qual_issn(self):
        """Features derived rom MeSH descriptor, MeSH qualifier, and Journal
        ISSN (abbreviated 'MQI')
        
        @return: Dictionary with 'mesh','qual','issn' keys, and with values
        being list of string feature of each type."""
        # Get MeSH headings, qualifiers and ISSN from article
        headings, quals = [], []
        for term in self.meshterms:
            headings.append(term[0])
            if len(term) > 1:
                for q in term[1:]:
                    if q not in quals:
                        quals.append(q)
        issns = [] if self.issn is None else [self.issn]
        # Get the feature vector while possibly them to the feature mapping
        return dict(mesh=headings, qual=quals, issn=issns)


    def feats_word(self):
        """Alphanumeric case-folded features derived from title and abstract/
        
        @return: Dictionary with the key 'word', and value being a list
        of word feature strings."""
        text = ""
        if self.title is not None:
            text = self.title + " "
        if self.abstract is not None:
            text += self.abstract
        # Collapse all non-alphabetics to single spaces and split into words
        words = re.sub(r'[^A-Za-z0-9-]+', ' ', text).split()
        # Keep only non-numeric, non-stop, length>1 words
        numbers = re.compile(r"^[0-9]+$")
        wordset = set(x for x in words if (len(x) > 1 and 
                x.lower() not in stopwords and not numbers.match(x)))
        return {"w":list(wordset)}
    
    def feats_word_folded(self):
        """Like feats_word, but with case-folding"""
        ft = self.feats_word()
        return {"w":list(set(x.lower() for x in ft["w"]))}
    
    def feats_word_nodash(self):
        """Like feats_word, but without hyphens"""
        text = ""
        if self.title is not None: text = self.title + " "
        if self.abstract is not None: text += self.abstract
        words = re.sub(r'[^A-Za-z0-9]+', ' ', text).split()
        numbers = re.compile(r"^[0-9]+$")
        wordset = set(x for x in words if (len(x) > 1 and x.lower() not in stopwords and not numbers.match(x)))
        return {"w":list(wordset)}
    
    def feats_word_num(self):
        """Like feats_word, but with numbers"""
        text = ""
        if self.title is not None: text = self.title + " "
        if self.abstract is not None: text += self.abstract
        words = re.sub(r'[^A-Za-z0-9-]+', ' ', text).split()
        wordset = set(x for x in words if (len(x) > 1 and x.lower() not in stopwords))
        return {"w":list(wordset)}


    def feats_author(self):
        """Space 'au', with features like 'JS Bloggs' for author names"""
        return {"a": [F+" "+L for F,L in self.authors if F and L]}


    def feats_wmqia(self):
        """Union of the L{feats_mesh_qual_issn}, L{feats_word} and
        L{feats_author} spaces. 
        @return: Feature dictionary with w,mesh,qual,issn,au keys."""
        features = self.feats_mesh_qual_issn()
        features.update(self.feats_word())
        features.update(self.feats_author())
        return features
    
    
    def feats_wmqia_filt(self):
        """Like L{feats_wmqia}, but remove word features that occur in MeSH
        features."""
        ft = self.feats_wmqia()
        meshtext = (" ".join(ft["mesh"] + ft["qual"])).lower()
        meshwords = set(re.sub(r'[^a-z0-9-]+', ' ', meshtext).split())
        ft["w"] = [x for x in ft["w"] if x.lower() not in meshwords]
        return ft


    def feats_iedb_word(self):
        """IEDB features from title and abstract only"""
        text = ""
        if self.title is not None:
            text = self.title + " "
        if self.abstract is not None:
            text += self.abstract
        words = set(x for x in re.sub(r'[^A-Za-z0-9-]+', ' ', text).split() 
                    if x.lower() not in stopwords)
        return {"iedb":list(words)}


    def feats_iedb_concat(self):
        """IEDB features: Concatenate title, abstract, MeSH, journal.  No 
        case folding."""
        text = []
        if self.title is not None:
            text.append(self.title)
        if self.abstract is not None:
            text.append(self.abstract)
        if self.journal is not None:
            text.append(self.journal)
        for first,last in self.authors:
            if last is not None:
                text.append(last)
        for term in self.meshterms:
            text.extend(term)
        text = " ".join(text)
        words = set(x for x in re.sub(r'[^A-Za-z0-9-]+', ' ', text).split()
                    if x.lower() not in stopwords)
        return {"iedb":list(words)}