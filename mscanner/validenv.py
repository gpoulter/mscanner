"""Environment for performing cross-validation-based analyses"""

                                     
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

from itertools import chain, izip
import logging as log
import numpy as nx
import os

from mscanner.configuration import rc
from mscanner.featuredb import FeatureDatabase
from mscanner.featuremap import FeatureMapping
from mscanner.plotting import Plotter
from mscanner import scorefile
from mscanner.scoring import FeatureInfo, countFeatures
from mscanner.support.gcheetah import TemplateMapper, FileTransaction
from mscanner.support.utils  import preserve_cwd, randomSample
from mscanner.validation import Validator, PerformanceStats

class ValidationEnvironment(object):
    """Manages cross-validation based analyses.
    
    @ivar featmap: Mapping feature ID <-> feature string, from L{__init__}
    @ivar featdb: Mapping doc ID -> list of feature IDs, from L{__init__}
    @ivar postfilter: Membership test for classifying things as positive, from L{__init__}

    @ivar positives: IDs of positive articles, from L{loadInputs}
    @ivar negatives: IDs of negative articles, from L{loadInputs}
    
    @ivar pscores: Scores of positive articles, from L{standardValidation}
    @ivar nscores: Scores of negative articles, from L{standardValidation}
    
    @ivar featinfo: FeatureInfo instance for calculating feature score, from L{writeReport}
    @ivar performance: Performance statistics based on article scores, from L{writeReport}
    """
    
    def __init__(self):
        """Constructor"""
        log.info("Loading article databases")
        self.featdb = FeatureDatabase(rc.featuredb, 'r')
        self.featmap = FeatureMapping(rc.featuremap)
        self.postFilterFunction = None
        
    @property
    def article_list(self):
        """Array with the PubMed IDs in the database 
        
        Array may be over 16 million members long, so the property and will
        take a while to load the first time it is accessed."""
        try:
            return self._article_list
        except AttributeError: 
            log.info("Loading article list")
            self._article_list = nx.array([int(x) for x in rc.articlelist.lines()])
            return self._article_list
        
    @preserve_cwd
    def standardValidation(self, pospath, negpath=None):
        """Loads data, performs validation, and writes report
        
        @keyword pospath: Location of input positive PMIDs
        
        @keyword negpath: Location of input negative PMIDs (or None,
        to select randomly from Medline)
        """
        # Construct report directory
        import time
        if rc.timestamp is None: rc.timestamp = time.time() 
        if not rc.valid_report_dir.exists():
            rc.valid_report_dir.makedirs()
            rc.valid_report_dir.chmod(0777)
        # FeatureInfo for this validation run
        self.featinfo = FeatureInfo(
            featmap = self.featmap, 
            pseudocount = rc.pseudocount, 
            mask = self.featmap.featureTypeMask(rc.exclude_types),
            frequency_method = rc.frequency_method,
            post_masker = rc.post_masker
        )
        os.chdir(rc.valid_report_dir)
        # Load saved results
        if rc.report_positives.isfile() and \
           rc.report_negatives.isfile():
            self.loadResults()
        # Calculate new results
        else:
            self.loadInputs(pospath, negpath)
            if len(self.positives) == 0:
                scorefile.no_valid_pmids_page(
                    list(scorefile.readPMIDs(rc.report_positives_broken)))
                return # nothing to do
            self.getResults()
            self.saveResults()
        self.writeReport()
        log.info("FINISHING VALIDATION %s", rc.dataset)
        rc.timestamp = None

    def loadInputs(self, pospath, negpath=None):
        """Load PubMed IDs from files.  

        @keyword pospath: Location of input positive PMIDs
        
        @keyword negpath: Location of input negative PMIDs (or None)"""
        log.info("Reading positive PubMed IDs")
        positives = list(scorefile.readPMIDs(pospath, include=self.featdb,
                        broken_name=rc.report_positives_broken))
        log.info("Reading negative PubMed IDs")
        if negpath is None:
            # No negatives provided
            # Clamp number of negatives if more are requested than are available
            maxnegs = len(self.article_list) - len(positives)
            if rc.numnegs > maxnegs:
                rc.numnegs = maxnegs
            # Take an appropriately sized sample of random citations
            self.negatives = randomSample(
                rc.numnegs, self.article_list, set(positives))
        else:
            # Read existing list of negatives from disk
            negatives = list(scorefile.readPMIDs(negpath, exclude=set(positives),
                         exclude_name=rc.report_negatives_exclude))
            self.negatives = nx.array(negatives, nx.int32)
        self.positives = nx.array(positives, nx.int32)
        
    @preserve_cwd
    def loadResults(self):
        """Load previously saved results instead of re-validating
        
        Assumes PubMed IDs in the files are decreasing by score, as written by
        L{scorefile.writePMIDScores}. """
        log.info("Loading result scores for %s", rc.dataset)
        os.chdir(rc.valid_report_dir)
        pscores, positives = zip(*scorefile.readPMIDs(rc.report_positives, withscores=True))
        nscores, negatives = zip(*scorefile.readPMIDs(rc.report_negatives, withscores=True))
        self.positives = nx.array(positives, nx.int32)
        self.pscores = nx.array(pscores, nx.float32)
        self.negatives = nx.array(negatives, nx.int32)
        self.nscores = nx.array(nscores, nx.float32)

    @preserve_cwd
    def saveResults(self):
        """Save validation scores to disk"""
        log.info("Saving result scores")
        os.chdir(rc.valid_report_dir)
        scorefile.writePMIDScores(rc.report_positives, izip(self.pscores, self.positives))
        scorefile.writePMIDScores(rc.report_negatives, izip(self.nscores, self.negatives))

    def getResults(self):
        """Calculate scores on citations using cross validation
        
        @note: All the self.featdb lookups are done beforehand, since lookups
        while busy with validation are slow."""
        s = self
        log.info("Performing cross-validation for for %s", rc.dataset)
        # Set up a postfilter (like geneDrugFilter) if necessary
        s.postfilter = None
        if s.postFilterFunction:
            s.postfilter = s.postFilterFunction()
        self.validator = Validator(
            featdb = dict((k,s.featdb[k]) for k in 
                          chain(s.positives,s.negatives)),
            featinfo = s.featinfo,
            positives = s.positives,
            negatives = s.negatives,
            nfolds = rc.nfolds,
            alpha = rc.alpha,
            postfilter = s.postfilter,
            )
        self.pscores, self.nscores = self.validator.validate()
        
    @preserve_cwd
    def geneDrugFilter(self):
        """Create a membership test for gene-drug association.
        
        To use the filter, assign this method to self.postFilterFunction
        
        @return: Set of PubMed IDs which have gene-drug co-occurrences.
        """
        from mscanner.pharmdemo.dbexport import (writeGeneDrugCountsCSV, 
                                                 countGeneDrug)
        from mscanner.pharmdemo.genedrug import getGeneDrugFilter
        log.info("Getting gene-drug associations") 
        os.chdir(rc.valid_report_dir)
        pos_arts = scorefile.getArticles(rc.articledb, self.positives)
        neg_arts = scorefile.getArticles(rc.articledb, self.negatives)
        gdfilter = getGeneDrugFilter(rc.genedrug, rc.drugtable, rc.gapscore)
        postfilter = set()
        for art in chain(pos_arts, neg_arts):
            gdresult = gdfilter(art)
            art.genedrug = gdresult
            if len(gdresult) > 0:
                postfilter.add(art.pmid)
        return postfilter
        #writeGeneDrugCountsCSV(countGeneDrug(pos_arts))

    @preserve_cwd
    def writeReport(self):
        """Write an HTML validation report. 
        
        @note: Only redraws figures for which output files do not already exist
        (likewise for term scores, but the index is always re-written).
        """
        os.chdir(rc.valid_report_dir)
        # Calculate feature score information
        self.featinfo.updateFeatureScores(
            pos_counts = countFeatures(
                len(self.featmap), self.featdb, self.positives),
            neg_counts = countFeatures(
                len(self.featmap), self.featdb, self.negatives),
            pdocs = len(self.positives),
            ndocs = len(self.negatives),
            )
        # Write term scores
        if not rc.report_term_scores.exists():
            log.debug("Writing term scores")
            import codecs
            f = codecs.open(rc.report_term_scores, "wb", "utf-8")
            try:
                self.featinfo.writeScoresCSV(f)
            finally:
                f.close()
        self.performance = PerformanceStats(
            self.pscores, self.nscores, rc.alpha)
        p = self.performance
        # Graph Plotting
        plotter = Plotter()
        overlap = None
        ##overlap, iX, iY = plotter.plotArticleScoreDensity(
        ##  rc.report_artscores_img, p.pscores, p.nscores, p.threshold)
        ##plotter.plotFeatureScoreDensity(
        ##  rc.report_featscores_img, elf.feature_info.scores)
        log.debug("Drawing graphs")
        # Article Score Histogram
        if not rc.report_artscores_img.exists():
            plotter.plotArticleScoreHistogram(
            rc.report_artscores_img, p.pscores, p.nscores, p.threshold)
        # Feature Score Histogram
        if not rc.report_featscores_img.exists():
            plotter.plotFeatureScoreHistogram(
            rc.report_featscores_img, self.featinfo.scores)
        # ROC Curve
        if not rc.report_roc_img.exists():
            plotter.plotROC(
            rc.report_roc_img, p.FPR, p.TPR, p.tuned.FPR)
        # Precision-Recall Curve
        if not rc.report_prcurve_img.exists():
            plotter.plotPrecisionRecall(
            rc.report_prcurve_img, p.TPR, p.PPV, p.tuned.TPR)
        # F-Measure Curve
        if not rc.report_fmeasure_img.exists():
            plotter.plotPrecisionRecallFmeasure(
            rc.report_fmeasure_img, p.uscores, p.TPR, p.PPV, 
            p.FM, p.FMa, p.threshold)
        # Index file
        log.debug("Writing index.html")
        mapper = TemplateMapper(root=rc.templates)
        ft = FileTransaction(rc.report_index, "w")
        tpl = mapper.validation(
            stats = self.featinfo.stats,
            linkpath = rc.templates.relpath().replace('\\','/') if rc.link_headers else None,
            overlap = overlap,
            p = self.performance,
            rc = rc,
            timestamp = rc.timestamp,
            notfound_pmids = list(scorefile.readPMIDs(rc.report_positives_broken)),
        ).respond(ft)
        ft.close()
