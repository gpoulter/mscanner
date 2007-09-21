"""Environment for performing query-based analyses."""

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

import codecs
import logging as log
import numpy as nx
from contextlib import closing

from mscanner.configuration import rc
from mscanner.support import dbshelve, utils
from mscanner import citationtable, cscore, featuredb, featuremap, scoring, scorefile


class Query:
    """Class for performing a single query

    @ivar env: L{Databases} to use.
    
    @ivar pmids: Set object with input PubMed IDs, from L{load_pmids}
    
    @ivar featinfo: FeatureInfo with feature scores, from L{query}
    
    @ivar inputs: List of (pmid, score) for input PMIDs, from L{make_results} or L{load_results}
    
    @ivar results: List of (pmid, score) for result PMIDs, from L{make_results} or L{load_results}
    """

    def __init__(self, outdir, env=None):
        """Initialise the query
        
        @param outdir: Path to directory for output files, which is created if
        it does not exist.
        
        @param env: L{Databases} to use.  If None we load them ourselves.
        """
        self.outdir = outdir
        if not outdir.exists():
            outdir.makedirs()
            outdir.chmod(0777)
        self.env = env if env else scorefile.Databases()
        self.pmids = None
        self.featinfo = None
        self.inputs = None
        self.results = None


    def load_pmids(self, input):
        """Generate the L{pmids} attribute from input, with validity checking
        
        @param input: List of PubMed IDs, or path to file containing the list.
        """
        if isinstance(input, set):
            self.pmids = input
        elif isinstance(input, basestring):
            log.info("Loading PubMed IDs from %s", input.basename())
            self.pmids = set(scorefile.read_pmids(
                input, include=self.env.featdb,
                broken_name=self.outdir/rc.report_input_broken))
            if len(self.pmids) == 0: # No valid PMIDs found
                scorefile.no_valid_pmids_page(
                    self.outdir/rc.report_index,
                    list(scorefile.read_pmids(self.outdir/rc.report_input_broken)))
        else:
            self.pmids = set(input)


    def make_featinfo(self):
        """Generate the L{featinfo} attribute using the L{pmids}
        as examples of relevant citations."""
        log.info("Calculating feature scores")
        pos_counts = scoring.count_features(
            len(self.env.featmap), self.env.featdb, self.pmids)
        self.featinfo = scoring.FeatureInfo(
            featmap = self.env.featmap,
            pos_counts = pos_counts,
            neg_counts = nx.array(self.env.featmap.counts, nx.int32) - pos_counts,
            pdocs = len(self.pmids),
            ndocs = self.env.featmap.numdocs - len(self.pmids),
            pseudocount = rc.pseudocount,
            mask = self.env.featmap.get_type_mask(rc.exclude_types),
            frequency_method = rc.frequency_method,
            post_masker = rc.post_masker)


    def query(self, input):
        """Performs a query given PubMed IDs as input
        
        @param input: Path to a list of PubMed IDs, or the list itself.
        """
        import time
        if rc.timestamp is None: 
            rc.timestamp = time.time() 
        self.load_pmids(input)
        if len(self.pmids) == 0: return
        self.make_featinfo()
        if (self.outdir/rc.report_input_scores).exists() \
           and (self.outdir/rc.report_result_scores).exists(): 
            self.load_results()
        else: 
            self.make_results()
            self.save_results()
        self.write_report()
        rc.timestamp = None # reset for next run
        log.info("FINISHING QUERY %s", rc.dataset)


    def load_results(self):
        """Read L{inputs} and L{results} from the report directory"""
        log.info("Loading saved results for %s", rc.dataset)
        self.inputs = list(scorefile.read_pmids(
            self.outdir/rc.report_input_scores, withscores=True))
        self.results = list(scorefile.read_pmids(
            self.outdir/rc.report_result_scores, withscores=True))


    def save_results(self):
        """Write L{inputs} and L{results} with scores in the report directory."""
        scorefile.write_scores(self.outdir/rc.report_input_scores, self.inputs)
        scorefile.write_scores(self.outdir/rc.report_result_scores, self.results)


    def make_results(self):
        """Perform the query to generate L{inputs} and L{results}"""
        log.info("Peforming query for dataset %s", rc.dataset)
        # Calculate score for each input PMID
        self.inputs = [ 
            (nx.sum(self.featinfo.scores[self.env.featdb[pmid]]), pmid)
            for pmid in self.pmids ]
        self.inputs.sort(reverse=True)
        # Calculate score for each result PMID
        self.results = list(cscore.score(
            rc.featurestream,
            self.env.featmap.numdocs,
            self.featinfo.scores,
            rc.limit, len(self.pmids),
            rc.threshold, self.pmids))
        self.results.sort(reverse=True)


    def write_report(self):
        """Write the HTML report for the query results
        
        @note: L{env.artdb} lookups are carried out beforehand because lookups
        while doing template output is extremely slow."""
        log.debug("Creating report for data set %s", rc.dataset)
        
        log.debug("Writing feature scores")
        with codecs.open(self.outdir/rc.report_term_scores, "wb", "utf-8") as f:
            self.featinfo.write_csv(f)

        log.debug("Writing input citations")
        self.inputs.sort(reverse=True)
        inputs = [ (s,self.env.artdb[str(p)]) for s,p in self.inputs]
        citationtable.writecitations(
            "input", inputs,
            self.outdir/rc.report_input_citations, 
            rc.citations_per_file)
        
        log.debug("Writing output citations")
        self.results.sort(reverse=True)
        outputs = [ (s,self.env.artdb[str(p)]) for s,p in self.results]
        citationtable.writecitations(
            "output", outputs, 
            self.outdir/rc.report_result_citations, 
            rc.citations_per_file)
        
        # Write ALL output citations to a single HTML, and a zip file
        outfname = self.outdir/rc.report_result_all
        zipfname = str(outfname + ".zip")
        citationtable.writecitations("output", outputs, outfname, len(outputs))
        from zipfile import ZipFile, ZIP_DEFLATED
        with closing(ZipFile(zipfname, "w", ZIP_DEFLATED)) as zf:
            zf.write(str(outfname), str(outfname.basename()))
        
        # Index.html
        log.debug("Writing index file")
        values = dict(
            stats = self.featinfo.stats,
            linkpath = rc.templates.relpath().replace('\\','/') if rc.link_headers else None,
            lowest_score = self.results[-1][0],
            num_results = len(self.results),
            best_tfidfs = self.featinfo.get_best_tfidfs(20),
            notfound_pmids = list(scorefile.read_pmids(self.outdir/rc.report_input_broken))
        )
        from Cheetah.Template import Template
        with utils.FileTransaction(self.outdir/rc.report_index, "w") as ft:
            Template(file=str(rc.templates/"results.tmpl"), 
                     filter="Filter", searchList=values).respond(ft)
