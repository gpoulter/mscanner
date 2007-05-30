"""Module implementing the Article object

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

class Article:
    """Stores an instance of a Medline citation.
    
    Instance variables are same as the constructor parameters. 
    
    @note: Article models a row
    (pmid,title,abstract,journal,issn,year,meshterms,authors), which 
    is stored in Berkeley DB to model a table keyed by PMID.

    @note: We avoided the OOP "God Object" problem where one class subsumes all
    behaviours which operate on its own instances in concert with other data
    source (such as Article.getScore(featurescores) and a million other
    methods).

    Passed via constructor:
        @ivar pmid: Integer PubMed ID or MEDLINE UI of the article.
        @ivar title: Title of the article
        @ivar abstract: Abstract of the article
        @ivar journal: Medline abbreviated journal title
        @ivar issn: ISSN code for the journal
        @ivar year: Year of publication
        @ivar meshterms: Set of Mesh terms associated with article.
        @ivar authors: Set of (initials,lastname) pairs of article authors
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
        """
        self.pmid = pmid
        self.title = title
        self.abstract = abstract
        self.journal = journal
        self.issn = issn
        self.year = year 
        self.meshterms = meshterms if meshterms else []
        self.authors = authors if authors else []

    def __repr__(self):
        import pprint
        pp = pprint.PrettyPrinter()
        s = "Article(pmid=%s,\n"+\
          "title=%s,\n"+\
          "abstract=%s,\n"+\
          "journal=%s\n"+\
          "issn=%s\n"+\
          "year=%s\n"+\
          "meshterms=%s\n"+\
          "authors=%s)\n"
        return s % (
            str(self.pmid),
            repr(self.title),
            repr(self.abstract),
            repr(self.journal),
            repr(self.issn),
            repr(self.year),
            pp.pformat(self.meshterms),
            pp.pformat(self.authors))
