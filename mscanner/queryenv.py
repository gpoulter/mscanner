"""Provides a high-level interface for performing query-based analyses.  

@note: This module relies heavily on rc parameters (many other modules like
scoring take these as function arguments instead)

                                   
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

import codecs
from itertools import chain
import logging as log
import numpy as nx
import os

from mscanner import dbshelve
from mscanner import statusfile
from mscanner.configuration import rc
from mscanner.featuredb import FeatureDatabase, FeatureStream
from mscanner.featuremap import FeatureMapping
from mscanner.gcheetah import TemplateMapper, FileTransaction
from mscanner.scoring import FeatureScoreInfo, filterDocuments
from mscanner.utils import (countFeatures, runMailer, readPMIDs, 
                            writePMIDScores, preserve_cwd)

class QueryEnvironment:
    """Class for managing MScanner analysis without passing around lots of
    parameters. RC parameters control the global data paths and output file
    names.

    @ivar featmap: Mapping between feature names and IDs
    @ivar featdb: Mapping from PubMed ID to list of features
    @ivar artdb: Mapping from PubMed ID to Article object
    @ivar narticles: Number of citations in PubMed
    
    @ivar input_pmids: List if PubMed IDs to query with
    @ivar feature_info: Object containing feature scores and statistics
    @ivar inputs: List of (pmid, score) for input PMIDs
    @ivar results: List of (pmid, score) for result PMIDs
    """
    
    def __init__(self):
        log.info("Loading article databases")
        self.featdb = FeatureDatabase(rc.featuredb, 'r')
        self.featmap = FeatureMapping(rc.featuremap)
        self.artdb = dbshelve.open(rc.articledb, 'r')
        self.narticles = int(rc.narticles.text())
        
    def __del__(self):
        log.info("Cleaning up")
        statusfile.close()
        self.featdb.close()
        self.artdb.close()
        
    def clearResults(self):
        del self.feature_info
        del self.inputs
        del self.results
        
    def loadInput(self, pmids_path):
        """Sets self.input_pmids"""
        log.info("Loading input PMIDs from %s", pmids_path)
        self.input_pmids = set(readPMIDs(pmids_path, include=self.featdb))
        statusfile.update(0, self.narticles - len(self.input_pmids))
        return self.input_pmids

    @preserve_cwd
    def standardQuery(self, pmids_path=None):
        """All-in-one for the most common operation

        @param pmids_path: Path to list of input PMIDs. If None, we assume
        loadInput() has already been called. """
        statusfile.start()
        if not rc.query_report_dir.exists():
            rc.query_report_dir.makedirs()
        if pmids_path is not None:
            self.loadInput(pmids_path)
        self.calculateFeatureScores()
        os.chdir(rc.query_report_dir)
        if rc.report_result_scores.exists() \
           and rc.report_positives.exists():
            self.loadSavedResults()
        else:
            self.performQuery()
        self.writeReport()
        statusfile.close()
        
    def calculateFeatureScores(self):
        """Calculate feature score information (self.feature_info)"""
        log.info("Calculating feature scores")
        pos_counts = countFeatures(
            len(self.featmap), self.featdb, self.input_pmids)
        neg_counts = nx.array(self.featmap.counts, nx.int32) - pos_counts
        self.feature_info = FeatureScoreInfo(
            pos_counts,  
            neg_counts,
            pdocs = len(self.input_pmids), 
            ndocs = self.narticles - len(self.input_pmids),
            pseudocount = rc.pseudocount, 
            featmap = self.featmap, 
            exclude_types = rc.exclude_types, 
            daniel = rc.dodaniel, 
            cutoff = rc.cutoff)
    
    @preserve_cwd
    def loadSavedResults(self):
        """Sets self.inputs and self.results by reading PMIDs from the report
        directory instead of actually querying."""
        os.chdir(rc.query_report_dir)
        log.info("Loading saved results for dataset %s", rc.dataset)
        self.inputs = list(readPMIDs(rc.report_positives, withscores=True))
        self.results = list(readPMIDs(rc.report_result_scores, withscores=True))
        
    @preserve_cwd
    def performQuery(self, save_results=True):
        """Perform the query.  If save_results, then list of pubmed IDs and 
        their scores are written in the report directory"""
        os.chdir(rc.query_report_dir)
        log.info("Peforming query for dataset %s", rc.dataset)
        # Calculate and write scores for each input PMID
        self.inputs = [ (pmid,nx.sum(self.feature_info.scores[self.featdb[pmid]])) 
                        for pmid in self.input_pmids ]
        if save_results:
            writePMIDScores(rc.report_positives, self.inputs)
        # Calculate and write score for each result PMID
        featstream = FeatureStream(file(rc.featurestream, "rb"))
        self.results = filterDocuments(
            featstream, self.feature_info.scores, 
            self.input_pmids, rc.limit, rc.threshold)
        if save_results:
            writePMIDScores(rc.report_result_scores, self.results)

    def gdFilter(self, export_pharmdemo=True):
        """Filter results and input articles for those containing gene-drug co-occurrences
        
        @param export_pharmdemo: If True, export results to PharmDemo database
        """
        from mscanner.genedrug import getGeneDrugFilter
        from mscanner.dbexport import exportDefault
        log.debug("Gene-drug associations on results")
        gdfilter = getGeneDrugFilter(rc.genedrug, rc.drugtable, rc.gapscore)
        self.gdarticles = []
        for pmid, score in chain(self.results, self.inputs):
            a = self.artdb[str(pmid)]
            a.genedrug = gdfilter(a)
            if len(a.genedrug) > 0:
                self.gdarticles.append(a)
        if export_pharmdemo:
            from mscanner.pharmdemo.dbexport import export_default
            export_default(self.gdarticles)
        return self.gdarticles

    @preserve_cwd
    def writeReport(self):
        """Write the HTML report for the query results
        """
        def writeCitations(mode, scores, fname, perfile):
            """Because 10000 citations in one HTML file is impossible to work
            with, function splits the results up into multiple files
            
            @param mode: 'input' or 'output'
    
            @param scores: List of (pmid, score) pairs in descending order of score
    
            @param fname: Base file name (e.g. results.html, resulting in
            results.html, results_02.html, ...)
    
            @param perfile: Number of citations per file (the very last file may
            have up to 2*perfile-1 citations)
            
            @note: Needed citations are first copied to an in-memory dict to
            work around a bug on mtree.stanford.edu where template takes
            forever to execute. This may be due to Cheetah doing something
            which takes forever on 22GB database. 
            """
            articles = dict()
            for pmid, score in scores:
                articles[str(pmid)] = self.artdb[str(pmid)]
            scores = list(scores)
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
                mapper.citations(
                    dataset = rc.dataset,
                    mode = mode, 
                    scores = towrite,
                    prev_file = fnames[count],
                    this_file = fnames[count+1],
                    next_file = fnames[count+2],
                    startfrom = start+1,
                    articles = articles).respond(
                        FileTransaction(fnames[count+1], "w"))
    
        log.debug("Writing report for data set %s", rc.dataset)
        mapper = TemplateMapper(root=rc.templates, kwargs=dict(filter="Filter"))
        if not rc.query_report_dir.isdir():
            rc.query_report_dir.mkdir()
        os.chdir(rc.query_report_dir)
        log.debug("Writing feature scores")
        self.feature_info.writeFeatureScoresCSV(codecs.open(rc.report_term_scores, "wb", "utf-8"))
        log.debug("Writing input citations")
        self.inputs.sort(key=lambda x:x[1], reverse=True)
        writeCitations("input", self.inputs, 
                       rc.report_input_citations, rc.citations_per_file)
        log.debug("Writing output citations")
        self.results.sort(key=lambda x:x[1], reverse=True)
        writeCitations("output", self.results, 
                       rc.report_result_citations, rc.citations_per_file)
        log.debug("Writing index file")
        index = mapper.results(
            rc = rc, 
            f = self.feature_info,
            num_results = len(self.results),
            lowest_score = self.results[-1][1])
        rc.report_index.write_text(str(index))
