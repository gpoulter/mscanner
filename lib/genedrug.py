#!/usr/bin/env python

"""Find associations of genes and drugs

@author: Graham Poulter
                                   

CachingGeneFinder -- Caches results of XMLRPC gene-finding queries
CachingGeneDrugLister -- Caches results of gene-drug association for articles
getGeneDrugFilter() -- Return a caching gene-drug association lister

"""

import re, xmlrpclib
import dbshelve
import logging as log
import os

class CachingGeneFinder:
    """Cache results of gene finding queries"""

    def __init__(self,cache):
        """Initialise cache
        
        Uses http://bionlp.stanford.edu/xmlrpc, function
        find_gene_and_protein_names, for the gene finder.

        @type cache: C{str:(str,int,int,float)}, or C{string}
        @param cache: Mapping of text to return values of
        L{findGenes}, or the file name of a persistent shelf
        containing such a mapping.
        """
        # Set up XMLRPC (possibly proxied)
        bionlp_uri="http://bionlp.stanford.edu/xmlrpc"
        proxy = os.environ.get("http_proxy", None)
        rpc_server = self.HTTPProxiedXMLRPC(bionlp_uri,proxy)
        self.geneFinder=rpc_server.find_gene_and_protein_names
        # Set up results cache
        if isinstance( cache, basestring ):
            self.cache = dbshelve.open( cache, 'c' )
        elif isinstance( cache, dict ):
            self.cache = cache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")

    def __del__(self):
        if isinstance( self.cache, dbshelve.Shelf):
            self.cache.close()

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
        key=repr(text)
        if key in self.cache:
            return self.cache[key]
        genes=self.geneFinder(text)
        result=[]
        for (name,start,end,score) in genes:
            if score >= 0.1: # throw away lowest scorers for efficiency
                result.append((name,start,end,score))
            self.cache[key] = result
        return result

    @staticmethod
    def HTTPProxiedXMLRPC(url,proxy):
        """Access an XMLRPC server from behind an HTTP proxy
        
        @type url: C{string}
        @param url: URL of XMLRPC service
        
        @type proxy: C{string}
        @param proxy: URL of HTTP proxy service
        """
        class ProxyTransport(xmlrpclib.Transport):
            def __init__(self,proxy):
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

class GeneDrugFilter:
    """Retrieve gene-drug associations for articles"""

    def __init__(self,drugs,geneFinder,threshold=0.85):
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
        self.fulldrugs=drugs
        self.drugs=self.stripDrugs(drugs)
        self.geneFinder=geneFinder
        self.threshold=threshold

    def __call__(self,article):
        """Automatically calls L{listGeneDrugsArticle}"""
        return self.listGeneDrugsArticle(article)

    @staticmethod
    def stripDrugs(drugs):
        """Return drug dictionary stripped down to lowercase characters

        @type drugs: C{{string:[string]}}
        @param drugs: Mapping from PharmGKB accession number (like
        PK12345) to drug aliases, where the first name in the list is
        preferred.

        @rtype: C{{string:[string]}}
        @return: Same as input except drug names consist of lowercase
        words separated by single space characters, all other
        characters having been stripped.
        """
        result = {}
        losepunct=re.compile(r'[0-9\s\!\@\#\$\%\^\&\*\(\){\}\[\]\-\:\;\"\'\?\/\.\,\<\>\`\~\-\_\=\+\\\|]+')
        for (PKID,names) in drugs.iteritems():
            lownames=[]
            for name in names:
                text=name.lower()
                text=losepunct.sub(' ',text).strip()
                lownames.append(text)
            result[PKID]=lownames
        return result

    def listDrugs(self,strings):
        """Return drugs found in each of a sequence of strings

        @type strings: C{[string]}
        @param strings: Strings to search for drugs.

        @rtype: C{[{string:string}]}
        @return: Mapping of PharmGKB ID to drug name for each string
        
        @note: Matching is approximate.  All non-alphabetics in both
        sentences and drugs are replaced with a single space before
        matching.
        """
        losepunct=re.compile(r'[0-9\s\!\@\#\$\%\^\&\*\(\){\}\[\]\-\:\;\"\'\?\/\.\,\<\>\`\~\-\_\=\+\\\|]+')
        strings=[" "+losepunct.sub(' ',s.lower())+" " for s in strings] # need extra spaces for later
        result=[dict() for s in strings]
        for (PKID,names) in self.drugs.iteritems():
            for (nidx,name) in enumerate(names):
                for (idx,s) in enumerate(strings):
                    if len(name)>4: # short names cause spurious matches
                        if s.find(" "+name+" ") != -1: # drug name must stand alone
                            result[idx][PKID]=self.fulldrugs[PKID][nidx]
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
        text=re.sub(r"\n",r" ",text)
        # Break on ellipsis
        text=re.sub(r"[.]{3}", r'&&&&&&&&&', text)
        # Break on exclamation, question
        text=re.sub(r'([!?]+\s+)', r'\1\n', text)
        # Break on full-stops followed by space, but not on
        # (title).(space) or .(space)(lowercase)
        titles=['Rev','Mr','Dr','Miss','Mrs','Ms','Sr','Prof','Ph']
        text=re.sub(r'('+'|'.join(titles)+r')[.]\s',r'\1&&& ',text)
        text=re.sub(r'[.](\s+[a-z])',r'&&&\1',text)
        text=re.sub(r'([.]\s+)',r'\1\n',text)
        # Split using newline markers, replace "&&&" with "."
        lines=[ re.sub(r'&&&',r'.',line) for line in text.split('\n') ]
        sentences=[]
        start_idx=0
        for line in lines:
            sentences.append((line,start_idx,start_idx+len(line)))
            start_idx+=len(line)
        return sentences

    def listGenes(self,text):
        """List genes found in text

        @type text: C{string}
        @param text: Input text to search for genes.

        @rtype: C{[(string,number,number)]}
        @return: (gene name,start index,end index) for each
        gene/protein name found.
        """
        genes=self.geneFinder(text)
        result=[]
        for (name,start,end,score) in genes:
            if score > self.threshold:
                result.append((name,start,end))
        return result

    def listGeneDrugs(self,text):
        """Return gene-drug co-occurrences in a text

        Search for genes in the text, and also split the text into
        sentences to search for drugs in each sentence.  Then
        aggregate co-occurrences of genes and drugs in each sentence.

        @type text: C{string}
        @param text: Text to parse, such as a title and abstract.

        @rtype: C{{string:[string]}}
        @return: Mapping from PharmGKB drug ID to genes co-occurring with
        the drug.
        """
        sentences=self.listSentences(text)
        genes_found=self.listGenes(text)
        drugs_found=self.listDrugs(s[0] for s in sentences)
        result={}
        for (idx,(sentence,sen_start,sen_end)) in enumerate(sentences):
            genes=set()
            # Add genes in this sentence
            for (gene,gene_start,gene_end) in genes_found:
                if gene_start >= sen_start and gene_end <= sen_end:
                    genes.add(gene)
            # Add drugs in this sentence
            if len(genes)!=0:
                for name in drugs_found[idx].itervalues():
                    if name not in result:
                        result[name]=genes
                    else:
                        result[name].update(genes)
        return result

    def listGeneDrugsArticle(self,article):
        """Return cached results from L{GeneDrugFilter.listGeneDrugs}

        @type article: C{Article}
        @param article: Article to cache result for.

        @rtype: C{string:set(string)}
        @return: The L{listGeneDrugs} result for an article's title
        and abstract joined together.
        """
        article.pmid=str(article.pmid)
        # DO NOT CHANGE NEXT FOUR LINES - WILL DISRUPT GENEFINDER CACHE
        text=article.title
        if text[-1:]!=".":
            text+=". "
        text+=article.abstract
        log.debug("Querying genes for article %s", article.pmid)
        return self.listGeneDrugs(text)

class CachingGeneDrugLister:

    def __init__( self, cache, gdfilter):
        """Caching proxy for gene-drug filtering
        @param cache: A mapping from PMIDs to a drug:[genes] mapping, or the filename for a shelf containing one. 
        @param gdfilter: Function which takes an article and returns a drug:[genes] mapping 
        """
        if isinstance(cache, basestring):
            self.cache = dbshelve.open(cache, 'c')
        elif isinstance(cache, dict):
            self.cache = cache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")
        self.gdfilter=gdfilter

    def __del__(self):
        if isinstance(self.cache, dbshelve.Shelf):
            self.cache.close()

    def __call__(self,article):
        """Automatically calls L{listGeneDrugs}"""
        return self.listGeneDrugs(article)

    def listGeneDrugs(self,article):
        """Return gene-drug association for an article"""
        key=str(article.pmid)
        if key in self.cache:
            return self.cache[key]
        result=self.gdfilter(article)
        self.cache[key]=result
        return result

def getGeneDrugFilter( gdcache, drugpickle, gapscorecache ):
    """Convenience function to create a fully caching gene-drug filter"""
    from cPickle import load
    return CachingGeneDrugLister(
        gdcache,
        GeneDrugFilter(
            load( file( drugpickle, "rb") ),
            CachingGeneFinder( gapscorecache )
        )
    )
