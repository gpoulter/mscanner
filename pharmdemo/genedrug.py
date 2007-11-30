"""Find co-occurrences of genes and drugs in abstract"""

from __future__ import with_statement
from __future__ import division

import logging
import os
import re
import xmlrpclib
from path import path

from mscanner.medline import Shelf


                                     
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


class DrugFinder:
    """Finds drug occurrences in text"""
    
    def __init__(self, source):
        """Initialise the drug finder
        
        @param source: Path to the drug table, or string with the contents.
        """
        if isinstance(source, path):
            with open(source, "rb") as f:
                source = f.read()
        self.drugtable = self._parse_drugs(source)
        self.patterns = self._make_patterns(self.drugtable)


    @staticmethod
    def _parse_drugs(text):
        """Parse table text into a drug dictionary
        
        @param text: String containing drug table to load.
        
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


    @staticmethod
    def _make_patterns(drugtable):
        """Create regular expressions for approximately matching any of the
        drug aliases in the drug dictionary.
        
        @param drugtable: Mapping from PharmGKB accession number (like PK12345) to
        list of drug aliases (the first name in the list is the preferred).
        
        @type drugtable: C{string:[string]}
        
        @return: Mapping from PharmGKB ID to a regex that matches any of the
        drug aliases (lowercased, with non-alphabetics stripped).
        
        @rtype: C{string:Pattern}
        """
        patterns = {}
        losechars = re.compile(r'[^a-z]+')
        for (PKID,names) in drugtable.iteritems():
            lownames = []
            for name in names:
                lname = losechars.sub(' ', name.lower()).strip()
                lname = lname.replace(' ',r'\s+')
                if lname not in lownames and len(lname) > 3:
                    lownames.append(lname)
            if len(lownames) > 0:
                patterns[PKID] = re.compile(r'\b(' + r'|'.join(lownames) + r')\b')
                #print patterns[PKID].pattern
        return patterns


    def find_drugs(self, text):
        """Return locations of drugs in text
        
        Matching is approximate. All non-alphabetics in both sentences and
        drugs are replaced with spaces before matching.
        
        @param text: Text to search for drugs.
        
        @return: List of (matched text, start index, end index, PharmGKB ID)
        
        @rtype: C{[(string,int,int,string)]}
        """
        # Remove punctuation, such that words are all space-separated
        reptext = re.sub(r'[^a-z0-9]'," ", text.lower())
        result = []
        for PKID, regex in self.patterns.iteritems():
            for m in regex.finditer(reptext):
                result.append((text[m.start():m.end()], m.start(), m.end(), PKID))
        return result



class ProxyTransport(xmlrpclib.Transport):
    """HTTP Proxy transport for use with xmlrpclib.ServerProxy
    
    @ivar proxy: string with hostname:port for HTTP proxy
    """
    def __init__(self, proxy, use_datetime=0):
        xmlrpclib.Transport.__init__(self, use_datetime)
        if proxy.startswith("http://"):
            proxy = proxy[len("http://"):]
        self.proxy = proxy

    def make_connection(self,host):
        self.realhost = host
        import httplib
        return httplib.HTTP(self.proxy)

    def send_request(self, connection, handler, request_body):
        connection.putrequest("POST", 'http://%s%s' % (self.realhost, handler))

    def send_host(self, connection, host):
        connection.putheader('Host', self.realhost)



class GeneFinder:
    """Perform gene-finding queries, with result caching.
    
    @ivar min_threshold: Minimum score for caching results
    
    @ivar threshold: Minimum score for returning results
    
    @ivar cache: Mapping from query text to GAPScore results.
    """

    def __init__(self, cache, threshold):
        """Initialise cache
        
        @note: Uses U{http://bionlp.stanford.edu/xmlrpc}, function
        C{find_gene_and_protein_names}.

        @param cache: Mapping from text to return values of L{gapscore}, or
        path to a persistent shelf to open."""
        self.min_threshold = 0.1
        self.threshold = threshold
        self._connect_bionlp()
        # Set up the results cache
        if isinstance(cache, path):
            self.cache = Shelf.open(cache, "c")
        elif isinstance(cache, dict):
            self.cache = cache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")


    def close(self):
        """Close the underlying cache if necessary"""
        if hasattr(self.cache, "close"):
            self.cache.close()


    def _connect_bionlp(self):
        """Connect to the BioNLP XMLRPC service (possibly proxied)"""
        bionlp_uri="http://bionlp.stanford.edu/xmlrpc"
        proxy = os.environ.get("http_proxy", None)
        if proxy is None:
            rpc_server = xmlrpclib.ServerProxy(bionlp_uri)
        else:
            rpc_server = xmlrpclib.ServerProxy(
                bionlp_uri, transport=ProxyTransport(proxy))
        self.find_gene_and_protein_names = rpc_server.find_gene_and_protein_names


    def gapscore(self, text):
        """Call the BioNLP GAPScore service with caching of results
        
        @param text: Text to search for genes and proteins
        
        @return: (gene, start, end, score) for each
        gene/protein name found.
        
        @rtype: C{[(str,int,int,float)]}
        """
        key = repr(text)
        try:
            return self.cache[key]
        except KeyError:
            genes = self.find_gene_and_protein_names(text)
            result = []
            for (name,start,end,score) in genes:
                if score >= self.min_threshold: 
                    result.append((name,start,end,score))
            self.cache[key] = result
            return result


    def find_genes(self, text):
        """List genes found in text
        
        @param text: Input text to search for genes.
        
        @return: (gene name, start index, end index) for each
        gene/protein name found.
        
        @rtype: C{[(string,number,number)]}
        """
        result = []
        for (name, start, end, score) in self.gapscore(text):
            if score > self.threshold:
                result.append((name, start, end))
        return result



class GeneDrugFinder:
    """Retrieve gene-drug associations for articles
    
    For caching results, use the object as a mapping from Article to
    drug:[genes]
    """

    def __init__(self, gdcache, genefinder, drugfinder):
        """Caching proxy for gene-drug filtering
        
        @param gdcache: A dictionary, or path to shelf. This maps PMID string
        to a drug:set(genes) mapping.
        
        @param genefinder: L{GeneFinder} instance
        
        @param drugfinder: L{DrugFinder} instance
        """
        self.genefinder = genefinder
        self.drugfinder = drugfinder
        if isinstance(gdcache, basestring):
            self.gdcache = Shelf.open(gdcache, 'c')
        elif isinstance(gdcache, dict):
            self.gdcache = gdcache
        else:
            raise ValueError("Cache is neither a filename nor dictionary")


    def close(self):
        """Close the gene-drug cache and gene finder"""
        if hasattr(self.gdcache, "close"):
            self.gdcache.close()
        self.genefinder.close()


    def __getitem__(self, article):
        """Get L{find_genedrugs_article} results from the cache"""
        key = str(article.pmid)
        try:
            result = self.gdcache[key]
        except KeyError:
            result = self.find_genedrugs_article(article)
        logging.debug("%s: %s", str(article.pmid), str(result))
        self.gdcache[key] = result
        return result


    @staticmethod
    def break_sentences(text):
        """Break text into a list of sentences
        @param text: Text to parse
        @return: (sentence, start index, end index) for sentences in text.
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


    def find_genedrugs_text(self, text):
        """Find gene-drug co-occurrences in a text
        
        Search for genes in the text, and split the text into sentences to
        search for drugs in each sentence. Then aggregate co-occurrences of
        genes and drugs in each sentence.
        
        @param text: Text to parse, such as a title and abstract.
        
        @return: Mapping from PharmGKB default drug name (obtained via
        PharmGKB ID) to genes co-occurring with the drug.
        
        @rtype: C{string:set([string])}
        """
        sentences = self.break_sentences(text)
        allgenes = self.genefinder.find_genes(text)
        alldrugs = self.drugfinder.find_drugs(text)
        # Filter out genes which overlap with a drug
        for drug, dstart, dfinish, pkid  in alldrugs:
            allgenes = [ (g,s,f) for g,s,f in allgenes if s >= dfinish or f < dstart ]
        result = {}
        for sentence, s_start, s_finish in sentences:
            genes = set([ g for g,s,f in allgenes if s >= s_start and f <= s_finish ])
            drugs = [ p for d,s,f,p in alldrugs if s >= s_start and f <= s_finish ]
            for pkid in drugs:
                if len(genes) > 0:
                    dname = self.drugfinder.drugtable[pkid][0].strip()
                    if dname not in result:
                        result[dname] = genes
                    else:
                        result[dname] |= genes
        return result


    def find_genedrugs_article(self, article):
        """Return gene-drug co-occurrences for an article
        
        @param article: Article object to find co-occurrences
        
        @return: Result of L{find_genedrugs_text} using concatenation of title
        and abstract."""
        article.pmid = str(article.pmid)
        # BEGIN DO_NOT_CHANGE
        text = article.title
        if text[-1:] != ".":
            text += ". "
        if article.abstract is not None:
            text += article.abstract
        # END DO_NOT_CHANGE
        result = self.find_genedrugs_text(text)
        return result



def open_genedrug_finder(genedrug_cache, drugtable, genefinder_cache):
    """Convenience function to create a fully caching gene-drug filter
    @param genedrug_cache: Dictionary, or path to DB of co-occurrences.
    @param drugtable: Path to table of drugs.
    @param genefinder_cache: Dictionary, or path to DB for gapscore results.
    """
    drugfinder = DrugFinder(drugtable)
    genefinder = GeneFinder(genefinder_cache, 0.85)
    return GeneDrugFinder(genedrug_cache, genefinder, drugfinder)
