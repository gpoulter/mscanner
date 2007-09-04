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
from mscanner import scoring, scorefile
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
        self.input_pmids = set(scorefile.readPMIDs(
            pmids_path, include=self.featdb,
            broken_name=rc.report_input_broken))
        return self.input_pmids
    
    def getFeatureInfo(self):
        """Calculate and return FeatureInfo object for current input"""
        log.info("Calculating feature scores")
        pos_counts = scoring.countFeatures(
            len(self.featmap), self.featdb, self.input_pmids)
        neg_counts = nx.array(self.featmap.counts, nx.int32) - pos_counts
        return scoring.FeatureInfo(
            featmap = self.featmap, 
            pos_counts = pos_counts, 
            neg_counts = neg_counts,
            pdocs = len(self.input_pmids), 
            ndocs = self.featmap.numdocs - len(self.input_pmids),
            pseudocount = rc.pseudocount, 
            mask = self.featmap.featureTypeMask(rc.exclude_types),
            frequency_method = rc.frequency_method,
            post_masker = rc.post_masker
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
        if pmids_path is not None: # Not already loaded PMIDs
            self.loadInput(pmids_path)
            if len(self.input_pmids) == 0: # No valid PMIDs found
                scorefile.no_valid_pmids_page(
                    list(scorefile.readPMIDs(rc.report_input_broken)))
                return
        # Perform query and write the report
        self.featinfo = self.getFeatureInfo()
        if rc.report_input_scores.exists() \
           and rc.report_result_scores.exists(): # Have results
            self.loadResults()
        else: # Generate new results
            self.getResults()
            self.saveResults()
        self.writeReport()
        rc.timestamp = None # reset for next run
        log.info("FINISHING QUERY %s", rc.dataset)
        
    def getGDFilterResults(self, pmids_path, export_db=False):
        """Filter results and input articles for those containing gene-drug
        co-occurrences
        
        @param export_pharmdemo: If True, export results to PharmDemo database
        """
        self.input_pmids = set(scorefile.readPMIDs(pmids_path, include=self.featdb))
        self.featinfo = self.getFeatureInfo()
        self.getResults()
        log.debug("Gene-drug associations on results")
        from mscanner.pharmdemo.genedrug import getGeneDrugFilter
        gdfilter = getGeneDrugFilter(rc.genedrug_db, rc.drugtable, rc.gapscore_db)
        gdarticles = []
        for score, pmid in chain(self.results, self.inputs):
            a = self.artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                gdarticles.append(a)
        gdfilter.close()
        if export_db == True:
            log.debug("Exporting database")
            from mscanner.pharmdemo.dbexport import GeneDrugExport
            gdexport = GeneDrugExport(gdarticles)
            gdexport.writeGeneDrugCountsCSV(rc.genedrug_csv)
            gdexport.exportText(rc.genedrug_sql)
        return gdarticles

    @preserve_cwd
    def loadResults(self):
        """Read self.inputs and self.results (PMIDs and scores) from
        the report directory"""
        os.chdir(rc.query_report_dir)
        log.info("Loading saved results for dataset %s", rc.dataset)
        self.inputs = list(scorefile.readPMIDs(rc.report_input_scores, withscores=True))
        self.results = list(scorefile.readPMIDs(rc.report_result_scores, withscores=True))
        
    @preserve_cwd
    def saveResults(self):
        """Write input/result PMIDs and scores to disk in report directory."""
        os.chdir(rc.query_report_dir)
        scorefile.writePMIDScores(rc.report_input_scores, self.inputs)
        scorefile.writePMIDScores(rc.report_result_scores, self.results)

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
            self.results = list(scoring.iterCScores(
                rc.cscore_path, rc.featurestream, 
                self.featmap.numdocs, self.featinfo.scores, 
                rc.limit, len(self.input_pmids),
                rc.threshold, self.input_pmids))
        else:
            # Read the feature stream using Python
            featurestream = FeatureStream(open(rc.featurestream, "rb"))
            self.results = scoring.iterScores(
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
        self.cumulative = scoring.retrievalTest(
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
        if not rc.query_report_dir.isdir():
            rc.query_report_dir.mkdir()
        os.chdir(rc.query_report_dir)
        # Feature Scores
        log.debug("Writing feature scores")
        import codecs
        f = codecs.open(rc.report_term_scores, "wb", "utf-8")
        self.featinfo.writeScoresCSV(f)
        f.close()
        # Write input citations
        from mscanner.citationtable import writecitations
        log.debug("Writing input citations")
        self.inputs.sort(reverse=True)
        inputs = [ (s,self.artdb[str(p)]) for s,p in self.inputs]
        writecitations("input", inputs,
                       rc.report_input_citations, rc.citations_per_file)
        # Write output citations
        log.debug("Writing output citations")
        self.results.sort(reverse=True)
        outputs = [ (s,self.artdb[str(p)]) for s,p in self.results]
        writecitations("output", outputs, 
                       rc.report_result_citations, rc.citations_per_file)
        # Write ALL output citations to a single HTML, and a zip file
        outfile = rc.report_result_all
        writecitations("output", outputs, outfile, len(outputs))
        from zipfile import ZipFile, ZIP_DEFLATED
        zf = ZipFile(str(outfile.basename() + ".zip"), "w", ZIP_DEFLATED)
        zf.write(str(outfile), str(outfile.basename()))
        zf.close()
        # Index.html
        log.debug("Writing index file")
        ft = FileTransaction(rc.report_index, "w")
        mapper = TemplateMapper(root=rc.templates, kwargs=dict(filter="Filter"))
        mapper.results(
            stats = self.featinfo.stats,
            linkpath = rc.templates.relpath().replace('\\','/') if rc.link_headers else None,
            lowest_score = self.results[-1][0],
            num_results = len(self.results),
            best_tfidfs = best_tfidfs,
            rc = rc, 
            timestamp = rc.timestamp,
            notfound_pmids = list(scorefile.readPMIDs(rc.report_input_broken)),
        ).respond(ft)
        ft.close()
