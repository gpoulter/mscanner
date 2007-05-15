"""Environment for constructing cross-validation-based analyses

ValidationEnvironment -- Environment for cross-validation

                                   
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
from itertools import chain, izip
import logging as log
import numpy as nx
import os

from mscanner import statusfile
from mscanner.configuration import rc
from mscanner.featuredb import FeatureDatabase
from mscanner.featuremap import FeatureMapping
from mscanner.gcheetah import TemplateMapper, FileTransaction
from mscanner.plotting import Plotter
from mscanner.scoring import FeatureInfo, countFeatures
from mscanner.utils import (getArticles, readPMIDs, runMailer, writePMIDScores,
                            preserve_cwd)
from mscanner.validation import Validator, PerformanceStats

class ValidationEnvironment:
    """Class to simplify constructing new cross-validation
    based analyses.
    
    From constructor:
    @ivar featmap: Mapping feature ID <-> feature string
    @ivar featdb: Mapping doc ID -> list of feature IDs
    @ivar featinfo: FeatureInfo instance for calculating feature score
    @ivar postfilter: Only classify positive if members of this

    From loadInputs():
    @ivar positives: IDs of positive articles
    @ivar negatives: IDs of negative articles
    
    From standardValidation():
    @ivar pscores: Scores of positive articles
    @ivar nscores: Scores of negative articles
    @ivar performance: Performance statistics based on article scores
    """
    
    def __init__(self):
        log.info("Loading article databases")
        self.featdb = FeatureDatabase(rc.featuredb, 'r')
        self.featmap = FeatureMapping(rc.featuremap)
        self.featinfo = FeatureInfo(
            self.featmap, 
            pseudocount = rc.pseudocount,
            mask = self.featmap.featureTypeMask(rc.exclude_types),
            daniel = rc.dodaniel, 
            cutoff = rc.cutoff)
        self.postfilter = None
        
    def __del__(self):
        log.debug("Cleaning up")
        statusfile.close()

    def reset(self):
        """Remove generated instance variables"""
        try:
            del self.positives
            del self.negatives
            del self.pscores
            del self.nscores
            del self.feature_info
            del self.performance
        except AttributeError:
            pass

    @preserve_cwd
    def standardValidation(self, pospath, negpath):
        """Loads data, performs validation, and writes report
        
        @param pospath: Location of input positive PMIDs
        @param negpath: Location of input negative PMIDs
        """
        # Report directory
        if not rc.valid_report_dir.exists():
            rc.valid_report_dir.makedirs()
        os.chdir(rc.valid_report_dir)
        statusfile.start(total=rc.nfolds)
        # Load saved results
        if rc.report_positives.isfile() and \
           rc.report_negatives.isfile():
            self.loadResults()
        # Calculate new results
        else:
            self.loadInputs(pospath, negpath)
            if rc.dogenedrug:
                self.doGeneDrug()
            self.getResults()
            self.saveResults()
        self.writeReport()
        statusfile.close()
        
    def loadInputs(self, pospath, negpath):
        """Load PubMed IDs from files"""
        log.info("Reading positive PubMed IDs")
        positives = list(readPMIDs(pospath, include=self.featdb))
        log.info("Reading negative PubMed IDs")
        negatives = list(readPMIDs(negpath, exclude=set(positives)))
        self.positives = nx.array(positives, nx.int32)
        self.negatives = nx.array(negatives, nx.int32)
        
    @preserve_cwd
    def loadResults(self):
        """Sets self.(positives, pscores, negatives, nscores),
        
        @note: writePMIDScores has written PMIDs decreasing by score
        """
        log.info("Loading result scores for %s", rc.dataset)
        os.chdir(rc.valid_report_dir)
        pscores, positives, = zip(*readPMIDs(rc.report_positives, withscores=True))
        nscores, negatives= zip(*readPMIDs(rc.report_negatives, withscores=True))
        self.positives = nx.array(positives, nx.int32)
        self.pscores = nx.array(pscores, nx.float32)
        self.negatives = nx.array(negatives, nx.int32)
        self.nscores = nx.array(nscores, nx.float32)

    @preserve_cwd
    def saveResults(self):
        """Save validation scores to disk"""
        os.chdir(rc.valid_report_dir)
        writePMIDScores(rc.report_positives, izip(self.pscores, self.positives))
        writePMIDScores(rc.report_negatives, izip(self.nscores, self.negatives))

    def getResults(self):
        """Calculate results by creating a Validator instance"""
        log.info("Performing cross-validation for for %s", rc.dataset)
        self.validator = Validator(
            featdb = self.featdb,
            featinfo = self.featinfo,
            positives = self.positives,
            negatives = self.negatives,
            nfolds = rc.nfolds,
            alpha = rc.alpha,
            postfilter = self.postfilter,
            )
        self.pscores, self.nscores = self.validator.validate()
        
    @preserve_cwd
    def doGeneDrug(self):
        from mscanner.pharmdemo.dbexport import writeGeneDrugCountsCSV, countGeneDrug
        from mscanner.pharmdemo.genedrug import getGeneDrugFilter
        log.info("Getting gene-drug associations") 
        os.chdir(rc.valid_report_dir)
        pos_arts = getArticles(rc.articledb, self.positives)
        neg_arts = getArticles(rc.articledb, self.negatives)
        gdfilter = getGeneDrugFilter(rc.genedrug, rc.drugtable, rc.gapscore)
        self.postfilter = set()
        for art in chain(pos_arts, neg_arts):
            gdresult = gdfilter(art)
            art.genedrug = gdresult
            if len(gdresult) > 0:
                self.postfilter.add(art.pmid)
        #writeGeneDrugCountsCSV(countGeneDrug(pos_arts))

    @preserve_cwd
    def writeReport(self):
        """Write an HTML validation report
        
        @note: sets self.performance and self.feature_info
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
            self.featinfo.writeScoresCSV(
                codecs.open(rc.report_term_scores, "wb", "utf-8"))
        self.performance = PerformanceStats(
            self.pscores, self.nscores, rc.alpha)
        p = self.performance
        # Write images
        plotter = Plotter()
        overlap = None
        ##overlap, iX, iY = plotter.plotArticleScoreDensity(
        ##  rc.report_artscores_img, p.pscores, p.nscores, p.threshold)
        ##plotter.plotFeatureScoreDensity(
        ##  rc.report_featscores_img, elf.feature_info.scores)
        if not rc.report_artscores_img.exists():
            plotter.plotArticleScoreHistogram(
            rc.report_artscores_img, p.pscores, p.nscores, p.threshold)
        if not rc.report_featscores_img.exists():
            plotter.plotFeatureScoreHistogram(
            rc.report_featscores_img, self.featinfo.scores)
        if not rc.report_roc_img.exists():
            plotter.plotROC(
            rc.report_roc_img, p.FPR, p.TPR, p.tuned.FPR)
        if not rc.report_prcurve_img.exists():
            plotter.plotPrecisionRecall(
            rc.report_prcurve_img, p.TPR, p.PPV, p.tuned.TPR)
        if not rc.report_fmeasure_img.exists():
            plotter.plotPrecisionRecallFmeasure(
            rc.report_fmeasure_img, p.pscores, p.TPR, p.PPV, 
            p.FM, p.FMa, p.threshold)
        # Write index file for validation output
        mapper = TemplateMapper(root=rc.templates)
        tpl = mapper.validation(
            stats = self.featinfo.getFeatureStats(),
            linkpath = rc.templates.relpath().replace('\\','/'),
            overlap = overlap,
            p = self.performance,
            rc = rc,
            timestamp = statusfile.timestamp,
        ).respond(FileTransaction(rc.report_index, "w"))
