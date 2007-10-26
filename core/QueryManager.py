"""Environment for performing query-based analyses."""

from __future__ import with_statement
from __future__ import division

import codecs
import logging as log
import numpy as nx
import time
from contextlib import closing
import warnings
warnings.simplefilter("ignore", UserWarning)

from mscanner.configuration import rc
from mscanner.medline import Shelf
from mscanner.medline.Databases import Databases
from mscanner.core.FeatureScores import FeatureScores, FeatureCounts
from mscanner.core import citationtable, cscore, iofuncs


                                     
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


class QueryManager:
    """Class for performing a single query

    @ivar env: L{Databases} to use.
    
    @ivar timestamp: Time at the start of the operation
    
    @ivar pmids: Set object with input PubMed IDs, from L{_load_pmids}
    
    @ivar featinfo: FeatureScores with feature scores, from L{query}
    
    @group From make_results or load_results: inputs, results
    
    @ivar inputs: List of (pmid, score) for input PMIDs
    
    @ivar results: List of (pmid, score) for result PMIDs
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
        self.timestamp = time.time()
        self.env = env if env else Databases()
        self.pmids = None
        self.featinfo = None
        self.inputs = None
        self.results = None
        # Use C score implementation if possible
        cscore.choose_score()


    def query(self, input):
        """Performs a query given PubMed IDs as input
        
        @param input: Path to a list of PubMed IDs, or the list itself.
        """
        try:
            self._load_input(input)
        except ValueError, e:
            log.error(str(e))
            iofuncs.no_valid_pmids_page(
                self.outdir/rc.report_index,
                list(iofuncs.read_pmids(self.outdir/rc.report_input_broken)))
            return
        self._make_feature_info()
        if (self.outdir/rc.report_input_scores).exists() \
           and (self.outdir/rc.report_result_scores).exists(): 
            self._load_results()
        else: 
            self._make_results()
            self._save_results()
        self._write_report()
        log.info("FINISHING QUERY %s", rc.dataset)


    def _load_input(self, input):
        """Generate the L{pmids} attribute from input, with validity checking
        
        @param input: List of PubMed IDs, or path to file containing the list.
        """
        if isinstance(input, set):
            self.pmids = input
        elif isinstance(input, basestring):
            log.info("Loading PubMed IDs from %s", input.basename())
            self.pmids = set(iofuncs.read_pmids(
                input, include=self.env.featdb,
                broken_name=self.outdir/rc.report_input_broken))
        else:
            self.pmids = set(input)


    def _make_feature_info(self):
        """Generate the L{featinfo} attribute using the L{pmids}
        as examples of relevant citations."""
        log.info("Calculating feature scores")
        self.featinfo = FeatureScores(
            featmap = self.env.featmap,
            pseudocount = rc.pseudocount,
            mask = self.env.featmap.get_type_mask(rc.exclude_types),
            make_scores = rc.make_scores,
            get_postmask = rc.get_postmask)
        pos_counts = FeatureCounts(
            len(self.env.featmap), self.env.featdb, self.pmids)
        self.featinfo.update(
            pos_counts = pos_counts,
            neg_counts = nx.array(self.env.featmap.counts, nx.int32) - pos_counts,
            pdocs = len(self.pmids),
            ndocs = self.env.featmap.numdocs - len(self.pmids))


    def _load_results(self):
        """Read L{inputs} and L{results} from the report directory"""
        log.info("Loading saved results for %s", rc.dataset)
        self.inputs = list(iofuncs.read_pmids(
            self.outdir/rc.report_input_scores, withscores=True))
        self.results = list(iofuncs.read_pmids(
            self.outdir/rc.report_result_scores, withscores=True))


    def _save_results(self):
        """Write L{inputs} and L{results} with scores in the report directory."""
        iofuncs.write_scores(self.outdir/rc.report_input_scores, self.inputs)
        iofuncs.write_scores(self.outdir/rc.report_result_scores, self.results)


    def _make_results(self):
        """Perform the query to generate L{inputs} and L{results}"""
        log.info("Peforming query for dataset %s", rc.dataset)
        # Calculate score for each input PMID
        self.inputs = zip(
            self.featinfo.scores_of(self.env.featdb, self.pmids), self.pmids)
        self.inputs.sort(reverse=True)
        # Calculate score for each result PMID
        self.results = list(cscore.score(
            rc.featurestream,
            self.env.featmap.numdocs,
            self.featinfo.scores,
            self.featinfo.offset,
            rc.limit, len(self.pmids),
            rc.threshold, self.pmids))
        self.results.sort(reverse=True)


    def _write_report(self):
        """Write the HTML report for the query results
        
        @note: Article database lookups are carried out beforehand because lookups
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
        if len(outputs) > 0:
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
            lowest_score = self.results[-1][0] if len(self.results) else 0,
            num_results = len(self.results),
            best_tfidfs = self.featinfo.get_best_tfidfs(20),
            timestamp = self.timestamp,
            notfound_pmids = list(iofuncs.read_pmids(self.outdir/rc.report_input_broken))
        )
        from Cheetah.Template import Template
        with iofuncs.FileTransaction(self.outdir/rc.report_index, "w") as ft:
            Template(file=str(rc.templates/"results.tmpl"), 
                     filter="Filter", searchList=values).respond(ft)