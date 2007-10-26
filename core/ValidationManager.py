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
from mscanner.core.Validator import CrossValidator


                                     
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



class SplitValidation:
    
    def __init__(self, outdir, env=None):
        """Constructor"""
        self.outdir = outdir
        if not outdir.exists():
            outdir.makedirs()
            outdir.chmod(0777)
        self.timestamp = time.time() 
        self.env = env if env else Databases()
        self.pscores, self.nscores = None, None


    def validation(self, fptrain, fntrain, fptest, fntest):
        """Carry out split-sample validation"""
        s = self
        s.ptrain, s.ntrain = None, None
        s.ptest, s.ntest = None, None
        s.notfound_pmids = []
        s.ptrain, broke, excl = iofuncs.read_pmids_careful(
            fptrain, s.env.featdb)
        s.ptest, broke, excl = iofuncs.read_pmids_careful(
            fptest, s.env.featdb)
        s.ntrain, broke, excl = iofuncs.read_pmids_careful(
            fntrain, s.env.featdb, set(s.ptrain))
        s.ntest, broke, excl = iofuncs.read_pmids_careful(
            fntest, s.env.featdb, set(s.ptest))
        if len(s.ptrain)>0 and len(s.ptest)>0 \
           and len(s.ntrain)>0 and len(s.ntest)>0:
            s._init_featinfo()
            s._calc_test_scores()
            s._calc_performance()
            s._write_report()


    def _init_featinfo(self):
        """Initialise L{featinfo} for use in validation"""
        self.featinfo = FeatureScores(
            featmap = self.env.featmap, 
            pseudocount = rc.pseudocount, 
            mask = self.env.featmap.get_type_mask(rc.exclude_types),
            make_scores = rc.make_scores,
            get_postmask = rc.get_postmask)


    def _calc_test_scores(self):
        """Use training sample for feature scores, then calc scores for testing sample."""
        s = self
        s.featinfo.update(
            FeatureCounts(len(s.featinfo), s.env.featdb, s.ptrain),
            FeatureCounts(len(s.featinfo), s.env.featdb, s.ntrain),
            len(s.ptrain), len(s.ntrain))
        s.pscores = nx.array(s.featinfo.scores_of(s.env.featdb, s.ptest))
        s.nscores = nx.array(s.featinfo.scores_of(s.env.featdb, s.ntest))


    def _calc_performance(self):
        """Calculate performance statistics"""
        self.performance = PerformanceStats(
            self.pscores, self.nscores, rc.alpha, rc.utility_r)
        self.perfrange = PerformanceRange(
            self.pscores, self.nscores, self.nfolds, self.performance.threshold)


    def _write_report(self):
        """Write an HTML validation report. 
        
        Only redraws figures for which output files do not already exist
        (likewise for term scores, but the index is always re-written)."""
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
            notfound_pmids = self.notfound_pmids,
        )
        from Cheetah.Template import Template
        with iofuncs.FileTransaction(self.outdir/rc.report_index, "w") as ft:
            Template(file=str(rc.templates/"validation.tmpl"), 
                     filter="Filter", searchList=values).respond(ft)


class ValidationManager(SplitValidation):
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
    @ivar notfound_pmids: List of input PMIDs not found in the database
    
    @group From constructor: featmap, featdb
    @ivar featmap: Mapping feature ID <-> feature string
    @ivar featdb: Mapping doc ID -> list of feature IDs

    @group From _make_results or _load_results: pscores, nscores
    @ivar pscores: Scores of positive articles
    @ivar nscores: Scores of negative articles
    """
    
    def __init__(self, outdir, env=None):
        """Constructor"""
        SplitValidation.__init__(self, outdir, env)
        self.nfolds = None
        self.featinfo = None


    def validation(self, pospath, negpath, nfolds):
        """Loads data, performs validation, and writes report
        
        @keyword pospath: Location of input positive PMIDs
        
        @keyword negpath: Location of input negative PMIDs (or None,
        to select randomly from Medline)
        """
        # Keep our own number of folds variables
        self.nfolds = nfolds
        self.notfound_pmids = []
        self._init_featinfo()
        # Load saved results
        try:
            self.positives, self.pscores = iofuncs.read_scores_array(
                self.outdir/rc.report_positives)
            self.negatives, self.nscores = iofuncs.read_scores_array(
                self.outdir/rc.report_negatives)
        # Perform the cross validation
        except IOError:
            if not self._load_input(pospath, negpath):
                return
            self._make_results()
            iofuncs.write_scores(self.outdir/rc.report_positives,
                                 izip(self.pscores, self.positives))
            iofuncs.write_scores(self.outdir/rc.report_negatives, 
                                 izip(self.nscores, self.negatives))
        # Report on the results
        self._calc_performance()
        self._general_feature_scores()
        self._write_report()
        log.info("FINISHING VALIDATION %s", rc.dataset)


    def _load_input(self, pospath, negpath):
        """Try to load positive and negative PubMed IDs
        
        @param pospath: Location of input positive PMIDs
        @param negpath: Location of input negative PMIDs
        """
        if isinstance(pospath, basestring):
            log.info("Loading positive PubMed IDs from %s", pospath.basename())
            self.positives, self.notfound_pmids, exclude = \
                iofuncs.read_pmids_careful(pospath, self.env.featdb)
        else:
            self.positives = nx.array(pospath)
        log.info("Loading negative PubMed IDs")
        if negpath is None:
            # Clamp number of negatives to the number available
            maxnegs = len(self.env.article_list) - len(self.positives)
            if rc.numnegs > maxnegs: rc.numnegs = maxnegs
            # Take a sample of random citations
            self.negatives = self.make_random_subset(
                rc.numnegs, self.env.article_list, set(self.positives))
        elif isinstance(negpath, basestring):
            # Read list of negative PMIDs from disk
            self.negatives, notfound, exclude = iofuncs.read_pmids_careful(
                    negpath, self.env.featdb, set(self.positives))
            self.notfound_pmids = list(self.notfound_pmids) + list(notfound)
            iofuncs.write_pmids(
                self.outdir/rc.report_negatives_exclude, exclude)
        else:
            self.negatives = nx.array(negpath)
        # Performing sanity checks
        iofuncs.write_pmids(
            self.outdir/rc.report_input_broken, self.notfound_pmids)
        if len(self.positives)>0 and len(self.negatives)>0:
            return True
        else:
            log.warning("No valid PubMed IDs in input... writing error page")
            iofuncs.no_valid_pmids_page(
                self.outdir/rc.report_index, self.notfound_pmids)
            return False
        return True


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


    def _make_results(self):
        """Calculate L{pscores} and L{nscores} using cross validation
        
        All the L{featdb} lookups are done beforehand, as lookups while busy
        with validation are slow."""
        log.info("Performing cross-validation for for %s", rc.dataset)
        self.validator = CrossValidator(
            featdb = dict((k,self.env.featdb[k]) for k in 
                          chain(self.positives,self.negatives)),
            featinfo = self.featinfo,
            positives = self.positives,
            negatives = self.negatives,
            nfolds = self.nfolds)
        self.pscores, self.nscores = self.validator.validate()


    def _general_feature_scores(self):
        """Calculate feature scores using all citations"""
        self.featinfo.update(
            pos_counts = FeatureCounts(
                len(self.env.featmap), self.env.featdb, self.positives),
            neg_counts = FeatureCounts(
                len(self.env.featmap), self.env.featdb, self.negatives),
            pdocs = len(self.positives),
            ndocs = len(self.negatives))



