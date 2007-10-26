"""Environment for performing cross-validation-based analyses"""

from __future__ import with_statement
from __future__ import division

from itertools import chain, izip
import logging as log
import numpy as nx
import time

import warnings
warnings.simplefilter("ignore", UserWarning)

from mscanner.configuration import rc
from mscanner.medline.Databases import Databases
from mscanner.core import iofuncs
from mscanner.core.FeatureScores import FeatureScores, FeatureCounts
from mscanner.core.PerformanceStats import PerformanceStats
from mscanner.core.PerformanceRange import PerformanceRange
from mscanner.core.Plotter import Plotter
from mscanner.core.Validator import Validator


                                     
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



class ValidationManager(object):
    """Manages the cross validation process.
    
    Jobs include loading input, saving and loading results files,
    and generating the validation report and figures.
    
    @ivar performance: L{PerformanceStats} instance
    @ivar perfrange: L{PerformanceRange} instance
    @ivar featinfo: L{FeatureScores} for calculating feature score
    @ivar timestamp: Time at the start of the operation

    @ivar nfolds: Number of validation folds
    @ivar positives: IDs of positive articles (from L{_load_positives})
    @ivar negatives: IDs of negative articles (from L{_load_negatives})
    
    @group From constructor: featmap, featdb
    @ivar featmap: Mapping feature ID <-> feature string
    @ivar featdb: Mapping doc ID -> list of feature IDs

    @group From _make_results or _load_results: pscores, nscores
    @ivar pscores: Scores of positive articles
    @ivar nscores: Scores of negative articles
    """
    
    def __init__(self, outdir, env=None):
        """Constructor"""
        self.outdir = outdir
        if not outdir.exists():
            outdir.makedirs()
            outdir.chmod(0777)
        self.timestamp = time.time() 
        self.env = env if env else Databases()
        self.nfolds = None
        self.positives = None
        self.pscores = None
        self.negatives = None
        self.nscores = None
        self.performance = None
        self.featinfo = None


    def validation(self, pospath, negpath, nfolds):
        """Loads data, performs validation, and writes report
        
        @keyword pospath: Location of input positive PMIDs
        
        @keyword negpath: Location of input negative PMIDs (or None,
        to select randomly from Medline)
        """
        # Keep our own number of folds variables
        self.nfolds = nfolds
        # Configure the FeatureScores object
        self.featinfo = FeatureScores(
            featmap = self.env.featmap, 
            pseudocount = rc.pseudocount, 
            mask = self.env.featmap.get_type_mask(rc.exclude_types),
            make_scores = rc.make_scores,
            get_postmask = rc.get_postmask
        )
        # Load saved results
        if (self.outdir/rc.report_positives).isfile() and \
           (self.outdir/rc.report_negatives).isfile():
            try:
                self._load_results()
            except ValueError, e:
                log.error(str(e))
                return
        # Calculate new results
        else:
            try:
                self._load_positives(pospath)
                self._load_negatives(negpath)
            except ValueError:
                # Error: unable to read any PMIDs
                log.error(str(e))
                iofuncs.no_valid_pmids_page(
                    self.outdir/rc.report_index,
                    list(iofuncs.read_pmids(self.outdir/rc.report_positives_broken)))
                return
            self._make_results()
            self._save_results()
        self.performance = PerformanceStats(
            self.pscores, self.nscores, rc.alpha)
        self.perfrange = PerformanceRange(
            self.pscores, self.nscores, self.nfolds, 
            self.performance.threshold)
        self._write_report()
        log.info("FINISHING VALIDATION %s", rc.dataset)


    def _load_positives(self, pospath):
        """Read positive PubMed IDs
        
        @param pospath: Location of input positive PMIDs (or a list of PMIDs
        directly.
        """
        log.info("Loading positive PubMed IDs")
        if isinstance(pospath, basestring):
            log.info("Loading PubMed IDs from %s", pospath.basename())
            self.positives = nx.array(list(iofuncs.read_pmids(
                pospath, include=self.env.featdb,
                broken_name=self.outdir/rc.report_positives_broken)))
        else:
            # See if we can convert to array
            self.positives = nx.array(pospath)


    @staticmethod
    def make_random_subset(k, pool, exclude):
        """Choose a random subset of k articles from pool
        
        This is better than the usual algorithm when the pool is large (say, 16
        million items), we don't mind if the order of pool gets scrambled, and
        we have to exclude certain items from being selected.
        
        @param k: Number of items to choose from pool
        @param pool: Array of items to choose from (will be scrambled!)
        @param exclude: Set of items that may not be chosen
        @return: A new array of the chosen items
        """
        from random import randint
        import numpy as nx
        n = len(pool)
        assert 0 <= k <= n
        for i in xrange(k):
            # Non-selected items are in 0 ... n-i-1
            # Selected items are n-i ... n
            dest = n-i-1
            choice = randint(0, dest) # 0 ... n-i-1 inclusive
            while pool[choice] in exclude:
                choice = randint(0, dest)
            # Move the chosen item to the end, where so it will be part of the
            # selected items in the next iteration. Note: this works using single
            # items - it but would break with slices due to their being views into
            # the vector.
            pool[dest], pool[choice] = pool[choice], pool[dest]
        # Phantom iteration: selected are n-k ... n
        return nx.array(pool[n-k:])


    def _load_negatives(self, negpath=None):
        """Load L{negatives} from file
        
        @note: Requires that L{positives} have already been loaded
        
        @note: If L{negpath} is None, we randomly select C{rc.numnegs} records
        from the list of citations in Medline.
    
        @param negpath: Location of input negative PMIDs 
        """
        log.info("Loading negative PubMed IDs")
        if negpath is None:
            # Clamp number of negatives if more are requested than are available
            maxnegs = len(self.env.article_list) - len(self.positives)
            if rc.numnegs > maxnegs: 
                rc.numnegs = maxnegs
            # Take an appropriately sized sample of random citations
            self.negatives = self.make_random_subset(
                rc.numnegs, self.env.article_list, set(self.positives))
        else:
            # Read existing list of negatives from disk
            negatives = list(
                iofuncs.read_pmids(
                    negpath, exclude=set(self.positives),
                    exclude_name=self.outdir/rc.report_negatives_exclude))
            self.negatives = nx.array(negatives, nx.int32)
        
        
    def _load_results(self):
        """Load L{positives}, L{pscores}, L{negatives}, L{nscores} using saved
        result files, instead of redoing the validation.
        
        This is useful when the style of the output page is updated. We assume
        PubMed IDs in the files are decreasing by score, as written by
        L{iofuncs.write_scores}. """
        log.info("Loading result scores for %s", rc.dataset)
        # Read positives scores
        pscores, positives = zip(*iofuncs.read_pmids(
            self.outdir/rc.report_positives, withscores=True))
        self.positives = nx.array(positives, nx.int32)
        self.pscores = nx.array(pscores, nx.float32)
        # Read negatives scores
        nscores, negatives = zip(*iofuncs.read_pmids(
            self.outdir/rc.report_negatives, withscores=True))
        self.negatives = nx.array(negatives, nx.int32)
        self.nscores = nx.array(nscores, nx.float32)


    def _save_results(self):
        """Save validation scores to disk"""
        log.info("Saving result scores")
        iofuncs.write_scores(self.outdir/rc.report_positives, 
                             izip(self.pscores, self.positives), sorted=False)
        iofuncs.write_scores(self.outdir/rc.report_negatives, 
                             izip(self.nscores, self.negatives), sorted=False)


    def _make_results(self):
        """Calculate L{pscores} and L{nscores} using cross validation
        
        All the L{featdb} lookups are done beforehand, as lookups while busy
        with validation are slow."""
        log.info("Performing cross-validation for for %s", rc.dataset)
        self.validator = Validator(
            featdb = dict((k,self.env.featdb[k]) for k in 
                          chain(self.positives,self.negatives)),
            featinfo = self.featinfo,
            positives = self.positives,
            negatives = self.negatives,
            nfolds = self.nfolds)
        self.pscores, self.nscores = self.validator.validate()


    def _use_global_feature_scores(self):
        """Calculate feature scores using all citations"""
        self.featinfo.update(
            pos_counts = FeatureCounts(
                len(self.env.featmap), self.env.featdb, self.positives),
            neg_counts = FeatureCounts(
                len(self.env.featmap), self.env.featdb, self.negatives),
            pdocs = len(self.positives),
            ndocs = len(self.negatives))


    def _write_report(self):
        """Write an HTML validation report. 
        
        Only redraws figures for which output files do not already exist
        (likewise for term scores, but the index is always re-written)."""
        self._use_global_feature_scores()
        # Write term scores
        if not (self.outdir/rc.report_term_scores).exists():
            log.debug("Writing term scores")
            import codecs
            with codecs.open(self.outdir/rc.report_term_scores, "wb", "utf-8") as f:
                self.featinfo.write_csv(f)
        # Draw graphs
        p = self.performance
        plotter = Plotter()
        log.debug("Drawing graphs")
        ##overlap, iX, iY = plotter.plot_score_density(
        ##  self.outdir/rc.report_artscores_img, p.pscores, p.nscores, p.threshold)
        # Article score histogram
        if not (self.outdir/rc.report_artscores_img).exists():
            plotter.plot_score_histogram(
            self.outdir/rc.report_artscores_img, p.pscores, p.nscores, p.threshold)
        # Feature score histogram
        if not (self.outdir/rc.report_featscores_img).exists():
            plotter.plot_feature_histogram(
            self.outdir/rc.report_featscores_img, self.featinfo.scores)
        # ROC curve
        if not (self.outdir/rc.report_roc_img).exists():
            plotter.plot_roc(
            self.outdir/rc.report_roc_img, p.FPR, p.TPR, p.tuned.FPR)
        # Precision-recall curve
        if not (self.outdir/rc.report_prcurve_img).exists():
            plotter.plot_precision(
            self.outdir/rc.report_prcurve_img, p.TPR, p.PPV, p.tuned.TPR)
        # F-Measure curve
        if not (self.outdir/rc.report_fmeasure_img).exists():
            plotter.plot_fmeasure(
            self.outdir/rc.report_fmeasure_img, p.uscores, p.TPR, p.PPV, 
            p.FM, p.FMa, p.threshold)
        # Write index file
        log.debug("Writing index.html")
        values = dict(
            stats = self.featinfo.stats,
            linkpath = rc.templates.relpath().replace('\\','/') if rc.link_headers else None,
            timestamp = self.timestamp,
            p = self.performance,
            perfrange = self.perfrange,
            notfound_pmids = list(iofuncs.read_pmids(self.outdir/rc.report_positives_broken)),
        )
        from Cheetah.Template import Template
        with iofuncs.FileTransaction(self.outdir/rc.report_index, "w") as ft:
            Template(file=str(rc.templates/"validation.tmpl"), 
                     filter="Filter", searchList=values).respond(ft)
