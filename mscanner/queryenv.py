"""Provides a high-level interface for performing query-based analyses.  

@note: This module relies heavily on rc parameters (many other modules like
scoring take these as function arguments instead)

"""

from __future__ import division, with_statement, absolute_import

                                               
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

from itertools import chain
import logging as log
import numpy as nx
import os

from mscanner.configuration import rc
from mscanner.featuredb import FeatureDatabase, FeatureStream
from mscanner.featuremap import FeatureMapping
from mscanner.scorefile import readPMIDs, writePMIDScores, emptyInputPage
from mscanner.scoring import (FeatureInfo, countFeatures, 
                              iterScores, iterCScores, retrievalTest)
from mscanner.support import dbshelve
from mscanner.support.gcheetah import TemplateMapper, FileTransaction
from mscanner.support.utils  import preserve_cwd

class QueryEnvironment:
    """Class for managing MScanner analysis without passing around lots of
    parameters. RC parameters control the global data paths and output file
    names.

    From constructor:
        @ivar featdb: Mapping from PubMed ID to list of features
        @ivar featmap: FeatureMapping instance between feature names and IDs
        @ivar artdb: Mapping from PubMed ID to Article object
    
    From loadInput():
        @ivar input_pmids: Set (not list!) of PubMed IDs as input to the query
    
    From standardQuery():
        @ivar featinfo: Object containing feature scores and statistics
    
    From loadResults() or getResults()
        @ivar inputs: List of (pmid, score) for input PMIDs
        @ivar results: List of (pmid, score) for result PMIDs
    """
    
    def __init__(self):
        """Initialise QueryEnvironment
        """
        log.info("Loading article databases")
        self.featdb = FeatureDatabase(rc.featuredb, 'r')
        self.featmap = FeatureMapping(rc.featuremap)
        self.artdb = dbshelve.open(rc.articledb, 'r')

    def __del__(self):
        self.featdb.close()
        self.artdb.close()
        
    def clearResults(self):
        """Get rid of old results, that we may query again"""
        try:
            del self.featinfo
            del self.input_pmids
            del self.results
        except AttributeError:
            pass

    def loadInput(self, pmids_path):
        """Read input PubMed IDs from pmids_path"""
        log.info("Loading input PMIDs from %s", pmids_path.basename())
        self.input_pmids = set(readPMIDs(
            pmids_path, outbase=rc.query_report_dir/"input", include=self.featdb))
        return self.input_pmids
    
    def getFeatureInfo(self):
        """Calculate and return FeatureInfo object for current input"""
        log.info("Calculating feature scores")
        pos_counts = countFeatures(
            len(self.featmap), self.featdb, self.input_pmids)
        neg_counts = nx.array(self.featmap.counts, nx.int32) - pos_counts
        return FeatureInfo(
            featmap = self.featmap, 
            pos_counts = pos_counts, 
            neg_counts = neg_counts,
            pdocs = len(self.input_pmids), 
            ndocs = self.featmap.numdocs - len(self.input_pmids),
            pseudocount = rc.pseudocount, 
            mask = self.featmap.featureTypeMask(rc.exclude_types),
            getFrequencies = rc.getFrequencies,
            getPostMask = rc.getPostMask
        )
    
    @preserve_cwd
    def standardQuery(self, pmids_path=None):
        """All-in-one for the most common operation.
        
        @note: re-calculates things which use RC parameters that may
        have changed since the last run.
        
        @param pmids_path: Path to list of input PMIDs. If None, we assume
        loadInput() has already been called. """
        import time
        if rc.timestamp is None: 
            rc.timestamp = time.time() 
        if not rc.query_report_dir.exists():
            rc.query_report_dir.makedirs()
            rc.query_report_dir.chmod(0777)
        os.chdir(rc.query_report_dir)
        if pmids_path is not None:
            self.loadInput(pmids_path)
            if len(self.input_pmids) == 0:
                # No valid PubMed IDs in the input, write error page
                emptyInputPage(list(
                    readPMIDs(rc.query_report_dir/"input.broken.txt")))
                return
        # Carry out the query and write the report
        self.featinfo = self.getFeatureInfo()
        if rc.report_input_scores.exists() \
           and rc.report_result_scores.exists():
            self.loadResults()
        else:
            self.getResults()
            self.saveResults()
        self.writeReport()
        rc.timestamp = None
        log.info("FINISHING QUERY %s", rc.dataset)
        
    @preserve_cwd
    def loadResults(self):
        """Read self.inputs and self.results (PMIDs and scores) from
        the report directory"""
        os.chdir(rc.query_report_dir)
        log.info("Loading saved results for dataset %s", rc.dataset)
        self.inputs = list(readPMIDs(rc.report_input_scores, withscores=True))
        self.results = list(readPMIDs(rc.report_result_scores, withscores=True))
        
    @preserve_cwd
    def saveResults(self):
        """Write input/result PMIDs and scores to disk in report directory."""
        os.chdir(rc.query_report_dir)
        writePMIDScores(rc.report_input_scores, self.inputs)
        writePMIDScores(rc.report_result_scores, self.results)

    def getResults(self):
        """Perform the query to generate input and result scores"""
        log.info("Peforming query for dataset %s", rc.dataset)
        # Calculate score for each input PMID
        self.inputs = [ 
            (nx.sum(self.featinfo.scores[self.featdb[pmid]]), pmid) 
            for pmid in self.input_pmids ]
        self.inputs.sort(reverse=True)
        # Calculate score for each result PMID
        if "use_cscores":
            # Read the feature stream using cscore
            self.results = list(iterCScores(
                rc.cscore_path, rc.featurestream, 
                self.featmap.numdocs, self.featinfo.scores, 
                rc.limit, len(self.input_pmids),
                rc.threshold, self.input_pmids))
        else:
            # Read the feature stream using Python
            featurestream = FeatureStream(open(rc.featurestream, "rb"))
            self.results = iterScores(
                featurestream, self.featinfo.scores, rc.limit,
                rc.threshold, self.input_pmids)
            featurestream.close()
        
    @preserve_cwd
    def testRetrieval(self, pmids_path):
        """Splits the input into training and testing, see how many
        testing PMIDs the query returns.
        
        @note: Writes files with the input PMIDs and scores, result PMIDs
        and scores, test pmids, and cumulative number of test PMIDs
         as a function of rank.
        
        @return: Array with cumulative test PMIDs as function of rank
        """
        if not rc.query_report_dir.exists():
            rc.query_report_dir.makedirs()
        os.chdir(rc.query_report_dir)
        # Split the input into train and test sections
        def split_input():
            input = list(self.input_pmids)
            import random
            random.shuffle(input)
            subset_size = int(rc.retrieval_test_prop * len(items))
            return set(items[:subset_size]), set(items[subset_size:])
        self.loadInput(pmids_path)
        self.input_pmids, self.input_test = split_input()
        # Get feature info and result scores
        self.featinfo = self.getFeatureInfo()
        rc.threshold = None
        self.getResults()
        self.saveResults()
        # Test the results against the test PMIDs
        self.cumulative = retrievalTest(
            [p for s,p in self.results], self.input_test)
        # Write the number of TP at each rank
        rc.report_retrieval_stats.write_lines(
            [str(x) for x in self.cumulative])
        # Write test PMIDs
        rc.report_retrieval_test_pmids.write_lines(
            [str(x) for x in self.input_test])
        # Graph TP vs FP (FP = rank-TP)
        from mscanner.plotting import Plotter
        plotter = Plotter()
        plotter.plotRetrievalGraph(
            rc.report_retrieval_graph, self.cumulative, len(self.input_test))
        return self.cumulative
    
    def getGDFilterResults(self, pmids_path, pharmdemo_db=None):
        """Filter results and input articles for those containing gene-drug
        co-occurrences
        
        @param export_pharmdemo: If True, export results to PharmDemo database
        """
        from mscanner.pharmdemo.genedrug import getGeneDrugFilter
        from mscanner.pharmdemo.dbexport import exportDefault
        self.loadInput(pmids_path)
        self.featinfo = self.getFeatureInfo()
        self.getResults()
        log.debug("Gene-drug associations on results")
        gdfilter = getGeneDrugFilter(rc.genedrug, rc.drugtable, rc.gapscore)
        self.gdarticles = []
        for pmid, score in chain(self.results, self.inputs):
            a = self.artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                self.gdarticles.append(a)
        if pharmdemo_db is not None:
            exportDefault(pharmdemo_db, self.gdarticles)
        return self.gdarticles

    @preserve_cwd
    def writeReport(self):
        """Write the HTML report for the query results
        
        @note: Writing of citation files is optimised by doing all of the
        self.artdb lookups beforehand.  Somehow, performing the lookups
        while busy with the template code is terribly slow. 
        """
        # Extract the top TF-IDF terms as (termid, (term,type), tfidf)
        from heapq import nlargest
        best_tfidfs = nlargest(
            20, enumerate(self.featinfo.tfidf), key=lambda x:x[1])
        best_tfidfs = [ 
            (t, tfidf, self.featmap[t], self.featinfo.scores[t], 
             self.featinfo.pos_counts[t], self.featinfo.neg_counts[t])
            for t, tfidf in best_tfidfs ]
        # Set up template and directory
        log.debug("Writing report for data set %s", rc.dataset)
        mapper = TemplateMapper(root=rc.templates, kwargs=dict(filter="Filter"))
        if not rc.query_report_dir.isdir():
            rc.query_report_dir.mkdir()
        os.chdir(rc.query_report_dir)
        # Feature Scores
        log.debug("Writing feature scores")
        import codecs
        f = codecs.open(rc.report_term_scores, "wb", "utf-8")
        self.featinfo.writeScoresCSV(f)
        f.close()
        # Write Input Citations
        log.debug("Writing input citations")
        self.inputs.sort(reverse=True)
        writeCitations(mapper.citations, "input", 
                       [ (s,self.artdb[str(p)]) for s,p in self.inputs],
                       rc.report_input_citations, rc.citations_per_file)
        # Write Output Citations
        log.debug("Writing output citations")
        self.results.sort(reverse=True)
        writeCitations(mapper.citations, "output",
                       [ (s,self.artdb[str(p)]) for s,p in self.results],
                       rc.report_result_citations, rc.citations_per_file)
        # Index.html
        log.debug("Writing index file")
        ft = FileTransaction(rc.report_index, "w")
        index = mapper.results(
            stats = self.featinfo.stats,
            #linkpath = rc.templates.relpath().replace('\\','/'),
            lowest_score = self.results[-1][0],
            num_results = len(self.results),
            best_tfidfs = best_tfidfs,
            rc = rc, 
            timestamp = rc.timestamp,
            notfound_pmids = list(readPMIDs(rc.query_report_dir/"input.broken.txt")),
        ).respond(ft)
        ft.close()

def writeCitations(template, mode, scores, fname, perfile):
    """Because 10000 citations in one HTML file is impossible to work with,
    writeCitations() takes the results and splits them across multiple files.
    
    @param template: Template class, initialised with parameters and
    its respond method is called with a FileTransaction.
    
    @param mode: 'input' or 'output'

    @param scores: List of (score, Article) in descending order of score

    @param fname: Base file name (e.g. results.html -->
    results.html, results_02.html, results_03.html etc.)

    @param perfile: Number of citations per file (the very last file 
    may however have up to 2*perfile-1 citations)
    """
    fnames = [None] + [fname] + [ 
        (fname.namebase + ("_%02d" % x) + fname.ext)
        for x in range(2,2+int(len(scores)/perfile)) ] + [None]
    starts = range(0, len(scores), perfile)
    if len(starts)>1 and (len(scores)-starts[-1]) < perfile:
        del starts[-1]
    for count, start in enumerate(starts):
        if count < len(starts)-1:
            towrite = scores[start:start+perfile]
        else:
            towrite = scores[start:]
            fnames[count+2] = None
        ft = FileTransaction(fnames[count+1], "w")
        template(
            cite_table = citationTable(start+1, towrite),
            dataset = rc.dataset,
            #linkpath = rc.templates.relpath().replace('\\','/'),
            mode = mode, 
            num_citations = len(towrite),
            prev_file = fnames[count],
            this_file = fnames[count+1],
            next_file = fnames[count+2],
            ).respond(ft)
        ft.close()

def citationTable(startfrom, scores):
    """Generate HTML table for citations using ElementTree
    
    @param startfrom: Rank of the first article in scores

    @param scores: List of (score, Article) in decreasing order of score
    
    @return: HTML for the <table> element containing citations
    
    @note: The philosophy here is to use Cheetah when there is more HTML than
    logic, and ElementTree when there is more logic than HTML. The old
    Cheetah template was getting cluttered from all the logic.  This 
    way also outputs less whitespace.
    
    @note: Original Cheetah template is a docstring at the end of the file.
    """
    from xml.etree.cElementTree import ElementTree, Element, SubElement
    table = Element("table", {"class":"cite"}, id="citations")
    ncols = 9 # Number of columns in the table
    cg = SubElement(table, "colgroup")
    SubElement(cg, "col", {"class":"classification"})
    SubElement(cg, "col", {"class":"rank"})
    SubElement(cg, "col", {"class":"score"})
    SubElement(cg, "col", {"class":"pmid"})
    SubElement(cg, "col", {"class":"year"})
    SubElement(cg, "col", {"class":"author"})
    SubElement(cg, "col", {"class":"abstract"})
    SubElement(cg, "col", {"class":"title"})
    SubElement(cg, "col", {"class":"journal"})
    thead = SubElement(table, "thead")
    tr = SubElement(thead, "tr")
    SubElement(tr, "th", title="Classification").text = "C"
    SubElement(tr, "th", title="Rank").text = "R"
    SubElement(tr, "th").text = "Score"
    SubElement(tr, "th").text = "PMID"
    SubElement(tr, "th").text = "Year"
    SubElement(tr, "th", title="Author").text = "Au"
    SubElement(tr, "th", title="Abstract").text = "Ab"
    SubElement(tr, "th").text = "Title"
    SubElement(tr, "th").text = "Journal"
    tbody = SubElement(table, "tbody")
    ncbi = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?"
    ncbi_pmid = ncbi+"cmd=Retrieve&db=pubmed&list_uids="
    ncbi_jour = ncbi+"CMD=search&DB=journals&term="
    for idx, (score, art) in enumerate(scores):
        pmid = str(art.pmid)
        tr = SubElement(tbody, "tr", {"class":"main"}, id=pmid)
        # Classification
        SubElement(tr, "td")
        # Rank
        SubElement(tr, "td").text = str(idx+startfrom)
        # Score
        SubElement(tr, "td").text = "%.2f" % score
        # PMID
        td = SubElement(tr, "td")
        SubElement(td, "a", href=ncbi_pmid+pmid).text = pmid
        # Year
        SubElement(tr, "td").text = str(art.year)
        # Expand Author button
        td = SubElement(tr, "td")
        td.text = "+" if art.authors else ""
        # Expand Abstract button
        td = SubElement(tr, "td")
        td.text = "+" if art.abstract else ""
        # Title
        td = SubElement(tr, "td", {"class":"title"})
        td.text = art.title.encode("utf-8")
        # ISSN
        td = SubElement(tr, "td")
        if art.issn:
            a = SubElement(td, "a", href=ncbi_jour+art.issn)
            a.text = art.journal if art.journal else art.issn
        # Expanded authors
        tr = SubElement(tbody, "tr", {"class":"author"}, id="au_"+pmid)
        if not art.authors:
            tr.set("style","display:none;")
        if art.authors:
            td = SubElement(tr, "td", colspan=str(ncols))
            td.text = ""
            for initials, lastname in art.authors:
                if initials: td.text += initials.encode("utf8") + " "
                if lastname: td.text += lastname.encode("utf8") + ", "
        # Expanded Abstract
        tr = SubElement(tbody, "tr", {"class":"abstract"}, id="ab_"+pmid)
        if not art.abstract:                        
            tr.set("style","display:none;")
        if art.abstract:
            td = SubElement(tr, "td", colspan=str(ncols))
            td.text = art.abstract.encode("utf8")
    from cStringIO import StringIO
    s = StringIO()
    ElementTree(table).write(s)
    return s.getvalue()
