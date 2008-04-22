"""Environment for performing cross-validation-based analyses"""

from __future__ import with_statement
from __future__ import division

import codecs
from contextlib import closing
from itertools import chain, izip
import logging
import numpy as nx
import random
import time

import warnings
warnings.simplefilter("ignore", UserWarning)

from mscanner.configuration import rc
from mscanner.medline.FeatureData import FeatureData
from mscanner.medline.FeatureVectors import random_subset
from mscanner.core import iofuncs
from mscanner.core.FeatureScores import FeatureScores
from mscanner.core.metrics import PerformanceVectors, PerformanceRange
from mscanner.core.Plotter import Plotter
from mscanner.core.Validator import cross_validate, count_features


                                     
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



class CrossValidation:
    """Carries out N-fold cross validation.
    
    @group Set in the constructor: outdir, dataset, fdata, timestamp, logfile

    @ivar outdir: Path to directory for output files, which is created if it
    does not exist.
    
    @ivar dataset: Title of the dataset to use when printing reports

    @ivar fdata: Selected FeatureData (document representation of articles).

    @ivar timestamp: Time at the start of the operation (set automatically by
    __init__).

    @ivar logfile: logging.FileHandler so that log output is
    sent to output directory as well (set automatically by __init__).


    @group Set by validation: featinfo, nfolds
    
    @ivar featinfo: L{FeatureScores} instance for calculating feature scores.

    @ivar nfolds: Number of cross validation folds (may not be relevant).


    @group Set by _load_input: positives, negatives, notfound_pmids, random_negatives
    
    @ivar positives: IDs of positive articles.

    @ivar negatives: IDs of negative articles.

    @ivar notfound_pmids: List of input PubMed IDs not found in the database.
    
    @group Set by _load_vectors: pos_vectors, neg_vectors, pos_counts, neg_counts
    pos_vectors: Feature vectors corresponding to positives.
    neg_vectors: Feature vectors corresponding to negatives.
    pos_counts: Number of occurrences of features in positives.
    neg_counts: Number of occurrences of features in negatives.
    
    @group Set by _crossvalid_scores: pscores, nscores
    @ivar pscores: Result scores for positive articles.
    @ivar nscores: Result scores for negative articles.

    @group Set by _get_performance: metric_vectors, metric_range
    @ivar metric_vectors: L{PerformanceVectors} instance
    @ivar metric_range: L{PerformanceRange} instance
    """

    
    def __init__(self, outdir, dataset, fdata):
        if not outdir.exists():
            outdir.makedirs()
            outdir.chmod(0777)
        # Attributes from constructor parameters
        self.outdir = outdir
        self.dataset = dataset
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
        self.random_negatives = False


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
        self.featinfo.numdocs = len(self.fdata.featuredb) # For scores_bgfreq
        # Try to load saved results
        try:
            logging.debug("Checking if there are saved results to load...")
            self.positives, self.pscores =\
                iofuncs.read_scores_array(self.outdir/rc.report_positives)
            self.negatives, self.nscores =\
                iofuncs.read_scores_array(self.outdir/rc.report_negatives)
            self._load_vectors(0)
            self.featinfo.update(self.pos_counts, self.neg_counts, 
                             len(self.positives), len(self.negatives))
            self.featinfo.stats
            logging.debug("Successfully loaded results")
        # Failed to load saved results, so perform cross validation
        except IOError:
            if not self._load_input(pos, neg):return
            self._load_vectors(rc.randseed)
            self._crossvalid_scores()
            iofuncs.write_scores(self.outdir/rc.report_positives,
                                 izip(self.pscores, self.positives))
            iofuncs.write_scores(self.outdir/rc.report_negatives, 
                                 izip(self.nscores, self.negatives))


    def report_validation(self):
        """Report cross validation results, using default threshold of 0"""
        if len(self.positives)>0 and len(self.negatives)>0:
            self._get_performance(0.0)
            self._write_report()

    
    def _load_vectors(self, randseed):
        """Shuffle the list of input articles (for the cross validator), then
        load corresponding feature vectors and derive counts for the PubMed IDs
        in L{positives} and L{negatives}.        
        @param randseed: Random seed for shuffling (0 for no shuffle, None for
        randomise)."""
        # Get the feature vectors for the positives (sorted by PubMed ID, then shuffle)
        pos_data = list((p,nx.array(v,nx.uint32)) for (p,d,v) in 
                        self.fdata.featuredb.get_records(self.positives))
        neg_data = list((p,nx.array(v,nx.uint32)) for (p,d,v) in 
                        self.fdata.featuredb.get_records(self.negatives))
        # Shuffle the feature vectors
        if randseed != 0:
            logging.debug("Shuffling pos/neg corpora with seed %s", str(randseed))
            random.seed(randseed)
            random.shuffle(pos_data)
            random.shuffle(neg_data)
        # Extract separate PubMed ID and vector lists
        self.positives, self.pos_vectors = zip(*pos_data)
        self.negatives, self.neg_vectors = zip(*neg_data)
        # Count feature occurrences
        self.pos_counts = count_features(len(self.fdata.featmap), self.pos_vectors)
        self.neg_counts = count_features(len(self.fdata.featmap), self.neg_vectors)


    def _crossvalid_scores(self):
        """Calculate article scores under cross validation. Use
        L{_load_vectors} to first shuffle the articles and load the required
        feature vectors. Before returning, we re-calculate feature scores using
        all of the training data.
        """
        logging.info("Calculating scores under cross validation.")
        self.pscores, self.nscores = cross_validate(
            self.featinfo, self.pos_vectors, self.neg_vectors, self.nfolds)
        logging.info("Updating feature scores after validation.")
        self.featinfo.update(self.pos_counts, self.neg_counts, 
                             len(self.positives), len(self.negatives))


    def _get_performance(self, threshold=None):
        """Calculate performance statistics.
        @param threshold: Score threshold, or None to use max F measure."""
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
        # Feature score histogram for selected features only
        plotter.plot_feature_histogram(
            self.outdir/rc.report_featscores_img, self.featinfo.scores[self.featinfo.selected])
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

        @param neg: Path to file of input negative PubMed IDs, or something
        convertible to integer array, or an integer representing the
        number of PubMed IDs to select at random from the database.
        
        @return: True if the load was successful, False otherwise.
        """
        if isinstance(pos, basestring):
            logging.info("Reading relevant PubMed IDs from %s", pos.basename())
            self.positives, self.notfound_pmids, exclude =\
                iofuncs.read_pmids_careful(pos, self.fdata.featuredb)
        else:
            logging.info("Using supplied array of positive PubMed IDs")
            self.positives = nx.array(pos, nx.int32)
        
        # Get random PubMed IDs
        if isinstance(neg, int):
            # Cap number of negatives requested to the maximum available
            maxnegs = len(self.fdata.featuredb) - len(self.positives)
            if neg > maxnegs: neg = maxnegs
            logging.info("Selecting %d random PubMed IDs for irrelevant examples." % neg)
            # Signal that the negatives are randomly sampled
            self.negatives = random_subset(
                neg, self.fdata.featuredb.pmids_array(), set(self.positives))
        # Read PubMed IDs from file
        elif isinstance(neg, basestring):
            logging.info("Reading irrelevant PubMed IDs from %s", neg.basename())
            self.negatives, notfound, exclude = iofuncs.read_pmids_careful(
                    neg, self.fdata.featuredb, set(self.positives))
            self.notfound_pmids = list(self.notfound_pmids) + list(notfound)
            iofuncs.write_pmids(self.outdir/rc.report_negatives_exclude, exclude)
        # Use PubMed IDs from parameter
        else:
            logging.info("Using supplied array of irrelevant PubMed IDs")
            self.negatives = nx.array(neg, nx.int32)
            
        # Report broken PubMed IDs
        iofuncs.write_pmids(self.outdir/rc.report_input_broken, self.notfound_pmids)
        
        # Checking that we have the input
        if len(self.positives)>0 and len(self.negatives)>0:
            return True
        else:
            logging.error("No valid PubMed IDs in at least one input (error page)")
            iofuncs.no_valid_pmids_page(self.outdir/rc.report_index, self.dataset, self.notfound_pmids)
            return False
        return True