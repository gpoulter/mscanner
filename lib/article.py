"""Module implementing the Article object

@author: Graham Poulter
                                        

Article -- Stores attributes of an article

"""

import pprint

class Article:
    """A simple wrapper for parsed Medline articles.
    
    Instance variables are same as the constructor parameters.

    """
    def __init__(self,
                 pmid=None,
                 title=None,
                 abstract=None,
                 journal=None,
                 issn=None,
                 year=None,
                 meshterms=None,
                 authors=None):
        """
        Initialise a new article with the given parameters.
        
        @param pmid: Integer PubMed ID or MEDLINE UI of the article.
        @param title: Title of the article
        @param abstract: Abstract of the article
        @param journal: Medline abbreviated journal title
        @param issn: ISSN code for the journal
        @param year: Year of publication
        @param meshterms: Set of Mesh terms associated with article.
        @param authors: Set of (initials,lastname) pairs of article authors
        """
        self.pmid = pmid
        self.title = title
        self.abstract = abstract
        self.journal = journal
        self.issn = issn
        self.year = year 
        self.meshterms = meshterms
        if meshterms is None:
            self.meshterms = list()
        self.authors = authors
        if authors is None:
            self.authors = list()

    def __repr__(self):
        pp = pprint.PrettyPrinter()
        astr = "Article(pmid=%d,\ntitle=%s,\nabstract=%s,\njournal=%s\nissn=%s\nyear=%s\nmeshterms=%s\nauthors=%s)\n"
        return astr % (
            self.pmid,
            repr(self.title),
            repr(self.abstract),
            repr(self.journal),
            repr(self.issn),
            repr(self.year),
            pp.pformat(self.meshterms),
            pp.pformat(self.authors))

