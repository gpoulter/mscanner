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
    
            @param scores: List of (PMID, score) in descending order of score
    
            @param fname: Base file name (results.html -->
            results.html, results_02.html, results_03.html etc.)
    
            @param perfile: Number of citations per file (the very last file 
            may however have up to 2*perfile-1 citations)
            
            @note: Needed citations are first copied to an in-memory dict,
            because trying to read the database while busy writing to
            HTML files really slows things down.
            """
            scores = list(scores)
            articles = dict((str(p),self.artdb[str(p)]) for p,s in scores)
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
                    linkpath = rc.templates.relpath().replace('\\','/'),
                    dataset = rc.dataset,
                    mode = mode, 
                    scores = towrite,
                    prev_file = fnames[count],
                    this_file = fnames[count+1],
                    next_file = fnames[count+2],
                    cite_table = citationTable(start+1, scores, articles)
                    ).respond(FileTransaction(fnames[count+1], "w"))
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
            f = self.feature_info,
            linkpath = rc.templates.relpath().replace('\\','/'),
            lowest_score = self.results[-1][1],
            num_results = len(self.results),
            rc = rc, 
            timestamp = statusfile.timestamp,
        ).respond(FileTransaction(rc.report_index, "w"))

def citationTable(startfrom, scores, articles):
    """Generate HTML table for citations using ElementTree
    
    @param startfrom: Rank of the first article in scores
    @param scores: List of (PMID, score) in decreasing order of score
    @param articles: Mapping from PMID to Article object
    
    @note: The philosophy here is to use Cheetah when there is more HTML than
    logic, and ElementTree when there is more logic than HTML. The old
    Cheetah template was getting cluttered from all the logic.  This 
    way also outputs less whitespace.
    """
    from xml.etree.cElementTree import ElementTree, Element, SubElement
    from cStringIO import StringIO
    table = Element("table", {"class":"cite"}, id="citations")
    ncols = 9 if articles else 4
    cg = SubElement(table, "colgroup")
    SubElement(cg, "col", {"class":"classification"})
    SubElement(cg, "col", {"class":"rank"})
    SubElement(cg, "col", {"class":"score"})
    SubElement(cg, "col", {"class":"pmid"})
    if articles:
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
    if articles:
        SubElement(tr, "th").text = "Year"
        SubElement(tr, "th", title="Author").text = "Au"
        SubElement(tr, "th", title="Abstract").text = "Ab"
        SubElement(tr, "th").text = "Title"
        SubElement(tr, "th").text = "Journal"
    tbody = SubElement(table, "tbody")
    ncbi = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?"
    ncbi_pmid = ncbi+"cmd=Retrieve&db=pubmed&list_uids="
    ncbi_jour = ncbi+"CMD=search&DB=journals&term="
    for idx, (pmid, score) in enumerate(scores):
        pmid = str(pmid)
        art = articles[pmid] if articles else None
        tr = SubElement(tbody, "tr", {"class":"main"}, id=pmid)
        # Classification
        SubElement(tr, "td", id="c_"+pmid, onclick="classify('%s')"%pmid)
        # Rank
        SubElement(tr, "td", id="c_"+pmid).text = str(idx+startfrom)
        # Score
        SubElement(tr, "td", id="s_"+pmid).text = "%.2f" % score
        # PMID
        SubElement(tr, "td", href=ncbi_pmid+pmid).text = pmid
        if art:
            # Year
            SubElement(tr, "td").text = str(art.year)
            # Expand Author
            td = SubElement(tr, "td", onclick="toggle('au_%s')"%pmid)
            td.text = "+" if art.authors else ""
            # Expand Abstract
            td = SubElement(tr, "td", onclick="toggle('ab_%s')"%pmid)
            td.text = "+" if art.abstract else ""
            # Title
            td = SubElement(tr, "td", {"class":"title"}, id="t_"+pmid)
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
    s = StringIO()
    ElementTree(table).write(s)
    return s.getvalue()
    
"""
## Old Cheetah Template for Historical Interest
<table id="citations" class="cite">
#set $ncols = 4
#if $articles
#set $ncols += 5
#end if
<colgroup>
  <col class="classification"></col>
  <col class="rank"></col>
  <col class="score"></col>
  <col class="pmid"></col>
  #if $articles
  <col class="year"></col>
  <col class="author"></col>
  <col class="abstract"></col>
  <col class="title"></col>
  <col class="journal"></col>
  #end if
</colgroup>
<thead>
  <tr>
    <th title="Classification">C</th>
    <th title="Rank">R</th>
    <th>Score</th>
    <th>PMID</th>
    #if $articles
    <th>Year</th>
    <th title="Author">Au</th>
    <th title="Abstract">Ab</th>
    <th>Title</th>
    <th>Journal</th>
    #end if
  </tr>
</thead>
<tbody>
#for idx, (pmid, score) in enumerate($scores)
#set $pmid = str(pmid)
#if $articles
#set $art = $articles[pmid]
#else 
#set $art = None
#end if
<tr id="#echo pmid#" class="main">
<td id="c_#echo pmid#" onclick="classify('#echo pmid#')"></td>
<td>#echo idx+$startfrom#</td>
<td id="s_#echo pmid#">#echo "%.2f" % score#</td>
<td>
<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;dopt=AbstractPlus&amp;list_uids=#echo pmid#">#echo pmid#</a>
</td>
#if art
<td>#echo art.year#</td>
<td#if art.authors# onclick="toggle('au_#echo pmid#')">+#else#>-#end if#</td>
<td#if art.abstract# onclick="toggle('ab_#echo pmid#')">+#else#>-#end if#</td>
<td id="t_#echo pmid#" class="title"> #echo art.title.encode("ascii","xmlcharrefreplace")# </td>
<td>
#if art.issn 
<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?CMD=search&amp;DB=journals&amp;term=#echo art.issn#">
#if art.journal##echo art.journal##else##echo art.issn##end if#
</a>
#end if
</td>
#end if
</tr>
#if art
<tr class="author" id="au_#echo pmid#"#if not art.authors# style="display:none;"#end if#>
#if art.authors
<td colspan="#echo ncols#">
#for initials, lastname in art.authors##slurp
#if initials##echo initials.encode("ascii","xmlcharrefreplace")##end if# #slurp
#if lastname##echo lastname.encode("ascii","xmlcharrefreplace")##end if#, #slurp
#end for#
##<script type="text/javascript">hide('au_#echo pmid#')</script>
</td>
#end if
</tr>
<tr class="abstract" id="ab_#echo pmid#"#if not art.abstract# style="display:none;"#end if#>
#if art.abstract
<td colspan="#echo ncols#" >
#echo art.abstract.encode("ascii","xmlcharrefreplace")
##<script type="text/javascript">hide('ab_#echo pmid#')</script>
</td>
#end if
</tr>
#end if
#end for
</tbody>
</table>
"""
