"""Environment for performing cross-validation-based analyses"""

from __future__ import with_statement
from __future__ import division

import codecs
from contextlib import closing
from itertools import chain, izip
import logging
import numpy as nx
import time

import warnings
warnings.simplefilter("ignore", UserWarning)

from mscanner.configuration import rc
from mscanner.medline.FeatureData import FeatureData
from mscanner.medline.ArticleData import ArticleData
from mscanner.core import iofuncs
from mscanner.core.FeatureScores import FeatureScores, FeatureCounts
from mscanner.core.metrics import (PerformanceVectors, PerformanceRange, 
                                   PredictedMetrics)
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



class CrossValidation(object):
    """Carries out N-fold cross validation.
    
    @group Set in the constructor: outdir, dataset, adata, fdata, timestamp, logfile

    @ivar outdir: Path to directory for output files, which is created if it
    does not exist.
    
    @ivar dataset: Title of the dataset to use when printing reports

    @ivar adata: Selected ArticleData (article database).
    
    @ivar fdata: Selected FeatureData (document representation of articles).

    @ivar timestamp: Time at the start of the operation (set automatically by
    __init__).

    @ivar logfile: logging.FileHandler so that log output is
    sent to output directory as well (set automatically by __init__).


    @group Set by L{validation}: featinfo, nfolds, positives, negatives,
    notfound_pmids, pscores, nscores
    
    @ivar featinfo: L{FeatureScores} instance for calculating feature scores.

    @ivar nfolds: Number of cross validation folds (may not be relevant).

    @ivar positives: IDs of positive articles.

    @ivar negatives: IDs of negative articles.

    @ivar notfound_pmids: List of input PMIDs not found in the database.

    @ivar pscores: Result scores for positive articles.

    @ivar nscores: Result scores for negative articles.


    @group Set by L{_get_performance}: metric_vectors, metric_range
    
    @ivar metric_vectors: L{PerformanceVectors} instance

    @ivar metric_range: L{PerformanceRange} instance
    """

    
    def __init__(self, outdir, dataset, adata, fdata):
        if not outdir.exists():
            outdir.makedirs()
            outdir.chmod(0777)
        # Attributes from constructor parameters
        self.outdir = outdir
        self.dataset = dataset
        self.adata = adata
        self.fdata = fdata
        self.timestamp = time.time() 
        self.logfile = iofuncs.open_logfile(self.outdir/rc.report_logfile)
        # Additional attributes for 
        self.nfolds = None
        self.pscores = None
        self.nscores = None
        self.notfound_pmids = []
        self.featinfo = None
        self.metric_vectors = None
        self.metric_range = None


    def __del__(self):
        iofuncs.close_logfile(self.logfile)


    def validation(self, pos, neg, nfolds=10):
        """Loads data and perform cross validation to calculate scores
        on that data.
        
        @note: This saves articles scores to the report directory, and 
        if possible it will load load those scores instead of calculating
        from scratch.
        
        @param pos, neg: Parameters for L{_load_input}
        
        @param nfolds: Number of validation folds to use.
        """
        logging.info("START: Cross validation for %s", self.dataset)
        # Keep our own number of folds attribute
        self.nfolds = nfolds
        self.notfound_pmids = []
        self.featinfo = FeatureScores.Defaults(self.fdata.featmap)
        self.featinfo.numdocs = self.adata.article_count # For scores_bgfreq
        # Try to load saved results
        try:
            self.positives, self.pscores = iofuncs.read_scores_array(self.outdir/rc.report_positives)
            self.negatives, self.nscores = iofuncs.read_scores_array(self.outdir/rc.report_negatives)
            self._update_featscores(self.positives, self.negatives)
            logging.debug("Successfully loaded saved cross validation results")
        # Failed to load saved results, so perform cross validation
        except IOError:
            if not self._load_input(pos, neg):
                return
            self.pscores, self.nscores = \
                self._crossvalid_scores(self.positives, self.negatives)
            iofuncs.write_scores(self.outdir/rc.report_positives,
                                 izip(self.pscores, self.positives))
            iofuncs.write_scores(self.outdir/rc.report_negatives, 
                                 izip(self.nscores, self.negatives))


    def report_validation(self):
        """Report cross validation results, using default threshold of 0"""
        if len(self.positives)>0 and len(self.negatives)>0:
            self._get_performance(0.0)
            self._write_report()

    
    def report_predicted(self, relevant_low, relevant_high, medline_size):
        """Experimental: report predicted query performance
        
        @param relevant_low: Minimum expected relevant articles in Medline
        
        @param relevant_high: Maximum expected relevant articles in Medline

        @param medline_size: Number of articles in rest of Medline, or None
        to use local database size minus relevant articles.
        """
        if len(self.positives)>0 and len(self.negatives)>0:
            logging.debug("Reporting performance prediction metrics")
            # Calculate the performance
            self._get_performance()
            if medline_size is None:
                medline_size = self.adata.article_count - len(self.positives)
            v = self.metric_vectors
            self.pred_low = PredictedMetrics(
                v.TPR, v.FPR, v.uscores, relevant_low, medline_size)
            self.pred_high = PredictedMetrics(
                v.TPR, v.FPR, v.uscores, relevant_high, medline_size)
            self._write_report()


    def _crossvalid_scores(self, positives, negatives):
        """Calculate article scores under cross validation.
        
        @note: L{positives} and L{negatives} come back shuffled by the cross
        validation partitioning! (use rc.randseed if you need the same shuffle
        every time)
        
        @note: Feature database lookups are slow so we cache them all
        in a dictionary (may use a lot of memory).
        
        @note: Before we return, we update L{featinfo} to use all of the
        training data (otherwise feature scores would reflect the
        last validation fold).
        
        @param positives: Vector of relevant PubMed IDs
        
        @param negatives: Vector of irrelevant PubMed IDs
        
        @return: Two vectors, containing scores for the positive and negative
        articles respectively, corresponding to the PubMed IDs in the
        (shuffled) L{positives} and L{negatives} vectors. 
        """
        logging.info("Calculating scores under cross validation.")
        self.validator = CrossValidator(
            featdb = dict((k,self.fdata.featuredb[k]) for k in 
                          chain(positives,negatives)),
            featinfo = self.featinfo,
            positives = positives,
            negatives = negatives,
            nfolds = self.nfolds,
            randseed = rc.randseed,
        )
        pscores, nscores = self.validator.validate()
        # Finally set feature scores using all available data
        self._update_featscores(positives, negatives)
        return pscores, nscores


    def _get_performance(self, threshold=None):
        """Calculate performance statistics.
        
        @param threshold: Specify a particular threshold, or None to estimate
        using F measure."""
        logging.info("Getting performance vectors (alpha=%s, utility_r=%s)", 
                     str(rc.alpha), str(rc.utility_r))
        v = PerformanceVectors(self.pscores, self.nscores, rc.alpha, rc.utility_r)
        self.metric_vectors = v
        if threshold is None:
            threshold, idx = v.threshold_maximising(v.FMa)
        else:
            threshold, idx = v.index_for(threshold)
        logging.info("Threshold is %f (uscores index=%d)", threshold, idx)
        average = v.metrics_for(idx)
        self.metric_range = PerformanceRange(
            self.pscores, self.nscores, self.nfolds, threshold, average)


    def _update_featscores(self, pos, neg):
        """Update the feature scores in L{featinfo} using the given 
        vectors of positive and negative citations."""
        logging.info("Updating feature scores.")
        self.featinfo.update(
            pos_counts = FeatureCounts(
                len(self.fdata.featmap), self.fdata.featuredb, pos),
            neg_counts = FeatureCounts(
                len(self.fdata.featmap), self.fdata.featuredb, neg),
            pdocs = len(pos),
            ndocs = len(neg))


    def _write_report(self):
        """Write an HTML validation report. 
        
        Only redraws figures for which output files do not already exist
        (likewise for term scores, but the index is always re-written)."""
        # Write term scores to file
        if not (self.outdir/rc.report_term_scores).exists():
            logging.debug("Writing features scores to %s", rc.report_term_scores)
            with closing(codecs.open(self.outdir/rc.report_term_scores, "wb", "utf-8")) as f:
                self.featinfo.write_csv(f, rc.max_output_features)
        # Aliases for the performance data
        p = self.metric_vectors
        t = self.metric_range.average
        # Do not overwriting existing plots
        plotter = Plotter(overwrite=False) 
        # Predicted precision/recall performance
        if hasattr(self, "pred_low") and hasattr(self, "pred_high"):
            plotter.plot_predictions(self.outdir/rc.report_prediction_img, 
                                     self.pred_low, self.pred_high)
        # Report cross validation results instead of prediction results
        else:
            # ROC curve
            plotter.plot_roc(
                self.outdir/rc.report_roc_img, p.FPR, p.TPR, t.FPR)
            # Precision-recall curve
            plotter.plot_precision(
                self.outdir/rc.report_prcurve_img, p.TPR, p.PPV, t.TPR)
            # F-Measure curve
            plotter.plot_fmeasure(
                self.outdir/rc.report_fmeasure_img, p.uscores, p.TPR, p.PPV, 
                p.FM, p.FMa, self.metric_range.threshold)
        # Article score histogram
        plotter.plot_score_histogram(
            self.outdir/rc.report_artscores_img, p.pscores, p.nscores, 
            self.metric_range.threshold)
        # Feature score histogram
        plotter.plot_feature_histogram(
            self.outdir/rc.report_featscores_img, self.featinfo.scores)
        # Write index file
        logging.debug("FINISH: Writing %s for %s", rc.report_index, self.dataset)
        from Cheetah.Template import Template
        with iofuncs.FileTransaction(self.outdir/rc.report_index, "w") as ft:
            Template(file=str(rc.templates/"validation.tmpl"), 
                     filter="Filter", searchList=dict(VM=self)).respond(ft)
   

    def _load_input(self, pos, neg):
        """Sets L{positives} and L{negatives} by various means
        
        @param pos: Path to file of input PubMed IDs, or something convertible
        to an integer array.

        @param neg: Path to file of input negative PMIDs, or something
        convertible to integer array, or an integer representing the
        number of PubMed IDs to select at random from the database.
        
        @return: True if the load was successful, False otherwise.
        """
        if isinstance(pos, basestring):
            logging.info("Reading positive PMIDs from %s", pos.basename())
            self.positives, self.notfound_pmids, exclude = \
                iofuncs.read_pmids_careful(pos, self.fdata.featuredb)
        else:
            logging.info("Using supplied array of positive PMIDs")
            self.positives = nx.array(pos, nx.int32)
            
        if isinstance(neg, int):
            logging.info("Selecting %d random PubMed IDs for background." % neg)
            maxnegs = self.adata.article_count - len(self.positives)
            if neg > maxnegs: neg = maxnegs
            self.negatives = self._random_subset(
                neg, self.adata.article_list, set(self.positives))
        elif isinstance(neg, basestring):
            logging.info("Reading negative PMIDs from %s", neg.basename())
            self.negatives, notfound, exclude = iofuncs.read_pmids_careful(
                    neg, self.fdata.featuredb, set(self.positives))
            self.notfound_pmids = list(self.notfound_pmids) + list(notfound)
            iofuncs.write_pmids(self.outdir/rc.report_negatives_exclude, exclude)
        else:
            logging.info("Using supplied array of negative PMIDs")
            self.negatives = nx.array(neg, nx.int32)
            
        # Writing out broken PubMed IDs
        iofuncs.write_pmids(
            self.outdir/rc.report_input_broken, self.notfound_pmids)
        
        # Checking that we have the input
        if len(self.positives)>0 and len(self.negatives)>0:
            return True
        else:
            logging.error("No valid PubMed IDs in at least one input (error page)")
            iofuncs.no_valid_pmids_page(
                self.outdir/rc.report_index, self.dataset, self.notfound_pmids)
            return False
        return True


    @staticmethod
    def _random_subset(k, pool, exclude):
        """Choose a random subset of k articles from pool
        
        This is a good algorithm when the pool is large (say, 16 million
        items), we don't mind if the order of pool gets scrambled, and we have
        to exclude certain items from being selected.
        
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
