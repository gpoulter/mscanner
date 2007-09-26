"""Find associations of genes and drugs

Members of this module:
    - L{parseDrugs}: Parse text version of PharmGKB drug table
    - L{GeneFinder}: Caches results of XMLRPC gene-finding queries
    - L{GeneDrugLister}: Use a gene-finder and drug table to find associations
    - L{GeneDrugListerCache}: Find associations for articles, caching results
    - L{getGeneDrugFilter}: Return a caching gene-drug association lister
"""

from __future__ import with_statement
from __future__ import division

                                     
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

import logging as log
import os
import re
import xmlrpclib

from mscanner.support import dbshelve


def parseDrugs(text):
    """Return a drug dictionary from specially formatted table

    @param text: String containing drug table to load.

    @see: Format of the drug table is described in README.html

    @raise ValueError: On incorrect formatting.

    @note: The specified file should map PharmGKB Accession ID's to drug names.  They are
    used to scan abstracts for occurrences of drugs.  One drug is
    described per line, and each line is to be formatted as follows::
        The drug's PharmGKB Accession ID
        [tab]
        The drug's main name in PharmGKB (might be comma-separated names)
        [tab]
        A list of the drug's generic names (if any), each terminated by '|'
        [tab]
        A list of the drug's trade names (if any), each terminated by '|'

    @return: Mapping from PharmGKB ID to list of aliases for the drug. The
    first alias is the preferred drug name."""
    drugs = {}
    text = text.replace("\r","")
    lines = text.split('\n')
    for line in lines:
        if line=="":
            continue
        # Extract tab-separated fields
        tabgroups = line.split("\t")
        if len(tabgroups) != 4:
            raise ValueError("Drug table: file must have four tab-separated fields. Bad line: <" + line + ">")
        PKID, main, trade, generic = tabgroups
        # Split trade names and generics
        if (trade != "" and trade[-1:] != "|") or \
           (generic != "" and generic[-1:] != "|"):
            raise ValueError("Drug table: generic/trade field did not end in '|'. Bad line: <"+line+">")
        # Split trade/generic lists using "|"
        trades = trade[:-1].split("|")
        generics = generic[:-1].split("|")
        drugs[PKID] = [ name for name in ([main]+trades+generics) if name != "" ]
    return drugs



class GeneFinder:
    """Perform gene-finding queries, with result caching.
    
    @ivar min_threshold: Minimum score for keeping results
    
    @ivar cache: Mapping from query texts to results"""

    def __init__(self, cache):
        """Initialise cache
        
        @note: Uses U{http://bionlp.stanford.edu/xmlrpc}, function
        C{find_gene_and_protein_names}.

        @param cache: Mapping from text to return values of L{findGenes}, or
        path to a persistent shelf to open."""
        # Default threshold
        self.min_threshold = 0.1
        # Set up XMLRPC (possibly proxied)
        bionlp_uri="http://bionlp.stanford.edu/xmlrpc"
        proxy = os.environ.get("http_proxy", None)
        if proxy is None:
            rpc_server = xmlrpclib.ServerProxy(bionlp_uri)
        else:
            if proxy.startswith("http://"):
                proxy = proxy[len("http://"):]
            rpc_server = self.HTTPProxiedXMLRPC(bionlp_uri, proxy)
        self.geneFinder=rpc_server.find_gene_and_protein_names
        # Set up results cache
        if isinstance(cache, basestring):
            self.cache = dbshelve.open(cache, "c")
        elif isinstance(cache, dict):
            self.cache = cache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")


    @staticmethod
    def HTTPProxiedXMLRPC(url, proxy):
        """Access an XMLRPC server from behind an HTTP proxy
        
        @param url: URL of XMLRPC service
        
        @param proxy: string with hostname:port for HTTP proxy
        """
        class ProxyTransport(xmlrpclib.Transport):
            def __init__(self, proxy, use_datetime=0):
                xmlrpclib.Transport.__init__(self, use_datetime)
                self.proxy = proxy
            def make_connection(self,host):
                self.realhost = host
                import httplib
                return httplib.HTTP(self.proxy)
            def send_request(self, connection, handler, request_body):
                connection.putrequest("POST", 'http://%s%s' % (self.realhost, handler))
            def send_host(self, connection, host):
                connection.putheader('Host', self.realhost)
        return xmlrpclib.ServerProxy(url,transport=ProxyTransport(proxy))


    def close(self):
        """Close the underlying cache"""
        if hasattr(self,"cache") and hasattr(self.cache, "close"):
            self.cache.close()
            del self.cache


    def __call__(self,text):
        """Shortcut for calling L{findGenes}"""
        return self.findGenes(text)


    def findGenes(self,text):
        """Finds genes in text, with caching of results

        @param text: Text to search for genes and proteins
        @type text: C{str}

        @return: (gene name,start index,end index,score) for each
        gene/protein name found.
        @rtype: C{[(str,int,int,float)]}
        """
        key = repr(text)
        if key in self.cache:
            return self.cache[key]
        genes = self.geneFinder(text)
        result = []
        for (name,start,end,score) in genes:
            if score >= self.min_threshold: 
                result.append((name,start,end,score))
            self.cache[key] = result
        return result



class GeneDrugLister:
    """Retrieve gene-drug associations for articles"""

    def __init__(self, drugs, geneFinder, threshold=0.85):
        """Initialise gene-drug filter

        @param geneFinder: Callable returning (gene name, start index, end
        index, score) for all genes/proteins located in the text. 
        
        @type geneFinder: C{f(string)->[(string,number,number,float)]}

        @param drugs: Mapping from PharmGKB ID's to drug names.  First
        name is used as the preferred name.

        @type drugs: C{string:[string]}
        
        @param threshold: Only use genes that score higher than this in
        geneFinder output.
        
        @type threshold: C{float}
        """
        self.fulldrugs = drugs
        self.drugs = self.stripDrugs(drugs)
        self.geneFinder = geneFinder
        self.threshold = threshold


    def __call__(self,article):
        """Shortcut for calling L{listGeneDrugsArticle}"""
        return self.listGeneDrugsArticle(article)


    def close(self):
        """Close the underlying gene-finding object"""
        if hasattr(self, "geneFinder") and hasattr(self.geneFinder, "close"):
            self.geneFinder.close()
            del self.geneFinder


    @staticmethod
    def stripDrugs(drugs):
        """Strip a drug dictionary down to lowercase characters and
        create a regular expression out of it.

        @param drugs: Mapping from PharmGKB accession number (like
        PK12345) to drug aliases (the first name in the list is the
        preferred).

        @type drugs: C{string:[string]}

        @return: Mapping from PharmGKB accession number to a regular expression
        to matches any of the drug names (lowercased, non-alphabetics
        stripped). 
        
        @rtype: C{string:Pattern} """
        result = {}
        losechars = re.compile(r'[^a-z]+')
        for (PKID,names) in drugs.iteritems():
            lownames = []
            for name in names:
                lname = losechars.sub(' ', name.lower()).strip()
                lname = lname.replace(' ',r'\s+')
                if lname not in lownames and len(lname) > 3:
                    lownames.append(lname)
            if len(lownames) > 0:
                result[PKID] = re.compile(r'\b(' + r'|'.join(lownames) + r')\b')
                #print result[PKID].pattern
        return result


    def listDrugs(self, text):
        """Return locations of drugs in text

        Matching is approximate. All non-alphabetics in both sentences and
        drugs are replaced with spaces before matching.
        
        @param text: Text to search for drugs.

        @return: List of (matched text, start index, end index, PharmGKB ID)
        
        @rtype: C{[(string,int,int,string)]}
        """
        # lose punctuation, with extra spaces so that every word boundary is marked by a space
        reptext = re.sub(r'[^a-z]'," ", text.lower())
        result = []
        for PKID, regex in self.drugs.iteritems():
            for m in regex.finditer(reptext):
                result.append((text[m.start():m.end()], m.start(), m.end(), PKID))
        return result


    @staticmethod
    def listSentences(text):
        """Return sentences from text
        
        @param text: Text to parse
        
        @return: (sentence, start index, end index) for sentences in text.
        
        @rtype: C{[(str,int,int)]}
        """
        # Remove existing newlines
        text = re.sub(r"\n",r" ",text)
        # Break on ellipsis
        text = re.sub(r"[.]{3}", r'&&&&&&&&&', text)
        # Break on exclamation, question
        text = re.sub(r'([!?]+\s+)', r'\1\n', text)
        # Break on full-stops followed by space, but not on
        # Mr.(space) or .(space)(lowercase)
        titles = ['Rev','Mr','Dr','Miss','Mrs','Ms','Sr','Prof','Ph']
        text = re.sub(r'('+'|'.join(titles)+r')[.]\s',r'\1&&& ',text)
        text = re.sub(r'[.](\s+[a-z])',r'&&&\1',text)
        text = re.sub(r'([.]\s+)',r'\1\n',text)
        # Split using newline markers, replace "&&&" with "."
        lines = [ re.sub(r'&&&',r'.',line) for line in text.split('\n') ]
        sentences = []
        start_idx = 0
        for line in lines:
            end_idx = start_idx+len(line)
            sentences.append((line, start_idx, end_idx))
            start_idx = end_idx
        return sentences


    def listGenes(self, text):
        """List genes found in text
        
        @param text: Input text to search for genes.
        
        @return: (gene name, start index, end index) for each
        gene/protein name found.
        
        @rtype: C{[(string,number,number)]}
        """
        genes = self.geneFinder(text)
        result = []
        for (name, start, end, score) in genes:
            if score > self.threshold:
                result.append((name, start, end))
        return result


    def listGeneDrugs(self, text):
        """Return gene-drug co-occurrences in a text
        
        Search for genes in the text, and also split the text into sentences to
        search for drugs in each sentence. Then aggregate co-occurrences of
        genes and drugs in each sentence.
        
        @param text: Text to parse, such as a title and abstract.
        
        @return: Mapping from PharmGKB default drug name (obtained via
        PharmGKB ID) to genes co-occurring with the drug.
        
        @rtype: C{string:set([string])}
        """
        sentences = self.listSentences(text)
        allgenes = self.listGenes(text)
        alldrugs = self.listDrugs(text)
        # Filter out genes which overlap with a drug
        for drug, dstart, dfinish, pkid  in alldrugs:
            allgenes = [ (g,s,f) for g,s,f in allgenes if s >= dfinish or f < dstart ]
        result = {}
        for sentence, s_start, s_finish in sentences:
            genes = set([ g for g,s,f in allgenes if s >= s_start and f <= s_finish ])
            drugs = [ p for d,s,f,p in alldrugs if s >= s_start and f <= s_finish ]
            for pkid in drugs:
                if len(genes) > 0:
                    dname = self.fulldrugs[pkid][0].strip()
                    if dname not in result:
                        result[dname] = genes
                    else:
                        result[dname] |= genes
        return result


    def listGeneDrugsArticle(self, article):
        """Get gene-drug co-occurrences for title/abstract of an article
        
        @param article: Article object to find co-occurrences
        
        @return: L{listGeneDrugs} result for concatenation of title
        and abstract.
        
        @rtype: C{string:set(string)}
        """
        article.pmid = str(article.pmid)
        # DO NOT CHANGE NEXT FIVE LINES (DISRUPTS GENEFINDER CACHE)
        text = article.title
        if text[-1:] != ".":
            text += ". "
        if article.abstract is not None:
            text += article.abstract
        result = self.listGeneDrugs(text)
        return result



class GeneDrugListerCache:

    def __init__(self, cache, gdfilter):
        """Caching proxy for gene-drug filtering
        
        @param cache: A mapping from PMIDs to a drug:[genes] mapping,
        or the path to a shelf for holding one.
        
        @param gdfilter: Function which takes an article and returns a
        drug:[genes] mapping.
        """
        if isinstance(cache, basestring):
            self.cache = dbshelve.open(cache, 'c')
        elif isinstance(cache, dict):
            self.cache = cache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")
        self.gdfilter=gdfilter


    def close(self):
        """Close the results cache and filter if possible"""
        if hasattr(self, "cache") and hasattr(self.cache, "close"):
            self.cache.close()
            del self.cache
        if hasattr(self, "gdfilter") and hasattr(self.gdfilter, "close"):
            self.gdfilter.close()
            del self.gdfilter


    def __call__(self, article):
        """Automatically calls L{listGeneDrugs}"""
        return self.listGeneDrugs(article)


    def listGeneDrugs(self, article):
        """Return gene-drug association for an article"""
        key = str(article.pmid)
        if key in self.cache:
            result = self.cache[key]
        else:
            result = self.gdfilter(article)
        log.debug("%s: %s", str(article.pmid), str(result))
        self.cache[key]=result
        return result



def getGeneDrugFilter(gdcache, drugtable, gapscorecache):
    """Convenience function to create a fully caching gene-drug filter
    
    @param gdcache: Dictionary or path to to DB for storing gene-drug co-occurrences
    
    @param drugtable: Path to table of drugs

    @param gapscorecache: Dictionary or path to DB for stroing gapscore outputs
    """
    with open(drugtable, "rb") as f:
        drugs = parseDrugs(f.read())
    genefinder = GeneFinder(gapscorecache)
    gdlister = GeneDrugLister(drugs, genefinder)
    return GeneDrugListerCache(gdcache, gdlister)
