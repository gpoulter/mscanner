"""Find associations of genes and drugs

parseDrugs() -- Parse text version of PharmGKB drug table
CachingGeneFinder -- Caches results of XMLRPC gene-finding queries
GeneDrugFilter -- Use a gene-finder and drug table to find associations
CachingGeneDrugLister -- Find associations for articles, caching results
getGeneDrugFilter() -- Return a caching gene-drug association lister

                                   

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

import logging as log
import os
import re
import unittest
import xmlrpclib

from mscanner import dbshelve

def parseDrugs(text):
    """Return a drug dictionary from specially formatted table

    @see: Format of the drug table is described in README.html

    @raise ValueError: On incorrect formatting.

    @type text: C{string}
    @param text: Drug table to load.

    @rtype: C{{string:[string]}}

    @note: The specified file should map PharmGKB Accession ID's to drug names.  They are
    used to scan abstracts for occurrences of drugs.  One drug is
    described per line, and each line is to be formatted as follows:
    
    The drug's PharmGKB Accession ID
    [tab]
    The drug's main name in PharmGKB (might be comma-separated names)
    [tab]
    A list of the drug's generic names (if any), each terminated by '|'
    [tab]
    A list of the drug's trade names (if any), each terminated by '|'

    @return: Mapping from PharmGKB ID to aliases for the drug.  The
    first alias is the preferred drug name.
    """
    drugs = {}
    lines = text.split('\n')
    mainsplit = re.compile(r"\A(\w+)\t([^\t]+)\t([^\t]*)\t(.*)\Z")
    subsplit = re.compile(r"\|")
    for line in lines:
        if line=="":
            continue
        # Check tab count in lines
        if line.count('\t') != 3:
            raise ValueError, "Drug lines must have three tabs (see README). Bad line was '"+line+"'"
        # Check splitting line into components
        (PKID,main,trade,generic) = mainsplit.match(line).groups()
        if trade != "" and trade[-1:] != "|":
            raise ValueError, "Drugs must terminate in '|' (see README). Bad line was '"+line+"'"
        if generic != "" and generic[-1:] != "|":
            raise ValueError, "Drugs must terminate in '|' (see README). Bad line was '"+line+"'"
        # Get names from each component and slurp into a list
        trades = subsplit.split(trade[:-1])
        generics = subsplit.split(generic[:-1])
        drugs[PKID] = [ name for name in [main]+trades+generics if name != "" ]
    return drugs


class CachingGeneFinder:
    """Cache results of gene finding queries
    
    @ivar min_threshold: Minimum score for keeping results
    
    @ivar cache: Mapping from query texts to results
    
    """

    def __init__(self, cache):
        """Initialise cache
        
        Uses http://bionlp.stanford.edu/xmlrpc, function
        find_gene_and_protein_names, for the gene finder.

        @type cache: C{str:(str,int,int,float)}, or C{string}
        @param cache: Mapping of text to return values of
        L{findGenes}, or the file name of a persistent shelf
        containing such a mapping.

        """
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

    def __del__(self):
        if hasattr(self,"cache") and isinstance(self.cache, dbshelve.Shelf):
            self.cache.close()
            del self.cache

    def __call__(self,text):
        """Automatically calls L{findGenes}"""
        return self.findGenes(text)

    def findGenes(self,text):
        """Finds genes in text, with caching of results

        @type text: C{str}
        @param text: Text to search for genes and proteins

        @rtype: C{[(str,int,int,float)]}
        @return: (gene name,start index,end index,score) for each
        gene/protein name found.
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

class GeneDrugFilter:
    """Retrieve gene-drug associations for articles"""

    def __init__(self, drugs, geneFinder, threshold=0.85):
        """Initialise gene-drug filter

        @type geneFinder: C{f(string)->[(string,number,number,float)]}
        @param geneFinder: Returns (gene name, start index, end index,
        score) for all genes/proteins located in the text.

        @type drugs: C{{string:[string]}}
        @param drugs: Mapping from PharmGKB ID's to drug names.  First
        name is used as the preferred name.

        @type threshold: C{float}
        @param threshold: Genes are only used if they score higher
        than this in the geneFinder output.
        """
        self.fulldrugs = drugs
        self.drugs = self.stripDrugs(drugs)
        self.geneFinder = geneFinder
        self.threshold = threshold

    def __call__(self,article):
        """Automatically calls L{listGeneDrugsArticle}"""
        return self.listGeneDrugsArticle(article)

    @staticmethod
    def stripDrugs(drugs):
        """Strip a drug dictionary down to lowercase characters and
        create a regular expression out of it.

        @type drugs: C{{string:[string]}}

        @param drugs: Mapping from PharmGKB accession number (like
        PK12345) to drug aliases (the first name in the list is the
        preferred).

        @rtype: C{{string:MatchObject}

        @return: Mapping from PharmGKB accession number to a regular
        expression which matches any of the drug names (lowercased,
        non-alphabetics stripped).
        """
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

        @param text: Text to search for drugs.

        @rtype: C{[(string,int,int,string)]}

        @return: List of (matched text, start index, end index, PharmGKB ID)
        
        @note: Matching is approximate.  All non-alphabetics in both
        sentences and drugs are replaced with spaces before matching.
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

        @type text: C{string}
        @param text: Text to parse

        @rtype: C{[(string,number,number)]}
        @return: (sentence,start index,end index) for sentences in the
        text.
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

        @rtype: C{[(string,number,number)]}
        
        @return: (gene name, start index, end index) for each
        gene/protein name found.
        """
        genes = self.geneFinder(text)
        result = []
        for (name, start, end, score) in genes:
            if score > self.threshold:
                result.append((name, start, end))
        return result

    def listGeneDrugs(self, text):
        """Return gene-drug co-occurrences in a text

        Search for genes in the text, and also split the text into
        sentences to search for drugs in each sentence.  Then
        aggregate co-occurrences of genes and drugs in each sentence.

        @param text: Text to parse, such as a title and abstract.

        @rtype: C{{string:set([string])}}

        @return: Mapping from PharmGKB default drug name (obtained via
        PharmGKB ID) to genes co-occurring with the drug.
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
        """Return cached results from L{GeneDrugFilter.listGeneDrugs}

        @param article: Article object to cache result for.

        @rtype: C{string:set(string)}
        
        @return: The L{listGeneDrugs} result for an article's title
        and abstract joined together.
        """
        article.pmid = str(article.pmid)
        # DO NOT CHANGE NEXT FOUR LINES (DISRUPTS GENEFINDER CACHE)
        text = article.title
        if text[-1:] != ".":
            text += ". "
        text += article.abstract
        log.debug("Querying genes for article %s", article.pmid)
        return self.listGeneDrugs(text)

class CachingGeneDrugLister:

    def __init__(self, cache, gdfilter):
        """Caching proxy for gene-drug filtering

        @param cache: A mapping from PMIDs to a drug:[genes] mapping,
        or the filename for a shelf containing one.

        @param gdfilter: Function which takes an article and returns a
        drug:[genes] mapping
        """
        if isinstance(cache, basestring):
            self.cache = dbshelve.open(cache, 'c')
        elif isinstance(cache, dict):
            self.cache = cache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")
        self.gdfilter=gdfilter

    def __del__(self):
        if hasattr(self, "cache") and isinstance(self.cache, dbshelve.Shelf):
            self.cache.close()
            del self.cache

    def __call__(self, article):
        """Automatically calls L{listGeneDrugs}"""
        return self.listGeneDrugs(article)

    def listGeneDrugs(self, article):
        """Return gene-drug association for an article"""
        key = str(article.pmid)
        if key in self.cache:
            return self.cache[key]
        result = self.gdfilter(article)
        log.debug("%s: %s", str(article.pmid), str(result))
        self.cache[key]=result
        return result

def getGeneDrugFilter(gdcache, drugtable, gapscorecache):
    """Convenience function to create a fully caching gene-drug filter"""
    return CachingGeneDrugLister(
        gdcache, GeneDrugFilter(
        parseDrugs(file(drugpickle, "rb").read()),
        CachingGeneFinder(gapscorecache)))
