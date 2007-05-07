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
from Gnuplot import Gnuplot
from itertools import chain, izip
import logging as log
import numpy as nx
import os

from mscanner import plotting
from mscanner import statusfile
from mscanner.configuration import rc
from mscanner.featuredb import FeatureDatabase
from mscanner.featuremap import FeatureMapping
from mscanner.gcheetah import TemplateMapper, FileTransaction
from mscanner.scoring import FeatureScoreInfo
from mscanner.utils import (countFeatures, getArticles, readPMIDs, 
                            runMailer, writePMIDScores, preserve_cwd)
from mscanner.validation import Validator, PerformanceStats

class ValidationEnvironment:
    """
    @ivar featmap: Mapping feature ID <-> feature string
    @ivar featdb: Mapping doc ID -> list of feature IDs
    @ivar narticles: Number of articles in Medline

    @ivar positives: IDs of positive articles
    @ivar negatives: IDs of negative articles
    @ivar pscores: Scores of positive articles
    @ivar nscores: Scores of negative articles
    @ivar feature_info: Feature scores using all positives & negatives
    @ivar performance: Performance statistics based on article scores
    """
    
    def __init__(self):
        log.info("Loading article databases")
        self.featdb = FeatureDatabase(rc.featuredb, 'r')
        self.featmap = FeatureMapping(rc.featuremap)
        self.narticles = int(rc.narticles.text())
        
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
        os.chdir(rc.valid_report_dir)
        statusfile.start(total=rc.nfolds)
        if not rc.valid_report_dir.exists():
            rc.valid_report_dir.makedirs()
        if rc.report_positives.isfile() and \
           rc.report_negatives.isfile():
            self.loadSavedScores()
        else:
            self.loadInputs(pospath, negpath)
            self.performValidation()
        self.writeReport()
        statusfile.close()
        
    def loadInputs(self, pospath, negpath):
        """Sets self.positives and self.negatives"""
        log.info("Reading positive PubMed IDs")
        positives = list(readPMIDs(pospath, include=self.featdb))
        log.info("Reading negative PubMed IDs")
        negatives = list(readPMIDs(negpath, exclude=set(positives)))
        self.positives = nx.array(positives, nx.int32)
        self.negatives = nx.array(negatives, nx.int32)
        return self.positives, self.negatives
        
    @preserve_cwd
    def loadSavedScores(self):
        """Sets self.(positives, pscores, negatives, nscores),
        
        @note: writePMIDScores put them on disk in decreasing order of score
        """
        log.info("Loading result scores for %s", rc.dataset)
        os.chdir(rc.valid_report_dir)
        positives, pscores = zip(*readPMIDs(rc.report_positives, withscores=True))
        negatives, nscores = zip(*readPMIDs(rc.report_negatives, withscores=True))
        self.positives = nx.array(positives, nx.int32)
        self.pscores = nx.array(pscores, nx.float32)
        self.negatives = nx.array(negatives, nx.int32)
        self.nscores = nx.array(nscores, nx.float32)
        return self.positives, self.pscores, self.negatives, self.nscores

    @preserve_cwd
    def performValidation(self):
        """Sets self.validator, self.pscores, self.nscores"""
        # Get which document ids have gene-drug assocations
        os.chdir(rc.valid_report_dir)
        log.info("Performing cross-validation for for %s", rc.dataset)
        self.genedrug_articles = None
        if rc.dogenedrug:
            self.doGeneDrug()
        self.validator = Validator(
            featmap = self.featmap,
            featdb = self.featdb,
            pos = self.positives,
            neg = self.negatives,
            nfold = rc.nfolds,
            pseudocount = rc.pseudocount,
            alpha = rc.alpha,
            daniel = rc.dodaniel,
            genedrug_articles = self.genedrug_articles,
            mask = self.featmap.featureTypeMask(rc.exclude_types)
            )
        self.pscores, self.nscores = self.validator.validate()
        writePMIDScores(rc.report_positives, izip(self.positives, self.pscores))
        writePMIDScores(rc.report_negatives, izip(self.negatives, self.nscores))
        
    @preserve_cwd
    def doGeneDrug(self):
        from mscanner.pharmdemo.dbexport import writeGeneDrugCountsCSV, countGeneDrug
        from mscanner.pharmdemo.genedrug import getGeneDrugFilter
        log.info("Getting gene-drug associations") 
        os.chdir(rc.valid_report_dir)
        self.genedrug_articles = set()
        pos_arts = getArticles(rc.articledb, self.positives)
        neg_arts = getArticles(rc.articledb, self.negatives)
        gdfilter = getGeneDrugFilter(rc.genedrug, rc.drugtable, rc.gapscore)
        for art in chain(pos_arts, neg_arts):
            gdresult = gdfilter(art)
            art.genedrug = gdresult
            if len(gdresult) > 0:
                genedrug_articles.add(art.pmid)
        writeGeneDrugCountsCSV(countGeneDrug(pos_arts))

    @preserve_cwd
    def writeReport(self):
        """Write an HTML validation report
        
        @note: sets self.performance and self.feature_info
        """
        os.chdir(rc.valid_report_dir)
        # Calculate feature score information
        self.feature_info = FeatureScoreInfo(
            pos_counts = countFeatures(
                len(self.featmap), self.featdb, self.positives),
            neg_counts = countFeatures(
                len(self.featmap), self.featdb, self.negatives),
            pdocs = len(self.positives),
            ndocs = len(self.negatives),
            pseudocount = rc.pseudocount,
            featmap = self.featmap,
            exclude_types = rc.exclude_types,
            daniel = rc.dodaniel,
            cutoff = rc.cutoff
            )
        if not rc.report_term_scores.exists():
            self.feature_info.writeFeatureScoresCSV(
                codecs.open(rc.report_term_scores, "wb", "utf-8"))
        overlap = None
        gp = Gnuplot(debug=1)
        if hasattr(self, "validator"):
            self.performance = self.validator.getPerformance()
        else:
            self.performance = PerformanceStats(
                self.pscores, self.nscores, rc.alpha)
        p = self.performance
        ##overlap, iX, iY = plotting.plotArticleScoreDensity(gp,
        ##  rc.report_artscores_img, p.pscores, p.nscores, p.tuned.threshold)
        ##plotting.plotFeatureScoreDensity(gp,
        ##  rc.report_featscores_img, elf.feature_info.scores)
        if not rc.report_artscores_img.exists():
            plotting.plotArticleScoreHistogram(gp, 
            rc.report_artscores_img, p.pscores, p.nscores, p.tuned.threshold)
        if not rc.report_featscores_img.exists():
            plotting.plotFeatureScoreHistogram(
            gp, rc.report_featscores_img, self.feature_info.scores)
        if not rc.report_roc_img.exists():
            plotting.plotROC(
            gp, rc.report_roc_img, p.FPR, p.TPR, p.tuned.FPR)
        if not rc.report_prcurve_img.exists():
            plotting.plotPrecisionRecall(
            gp, rc.report_prcurve_img, p.TPR, p.PPV, p.tuned.TPR)
        if not rc.report_fmeasure_img.exists():
            plotting.plotPrecisionRecallFmeasure(
            gp, rc.report_fmeasure_img, p.pscores, p.TPR, p.PPV, 
            p.FM, p.FMa, p.tuned.threshold)
        #Write index file for validation output
        mapper = TemplateMapper(root=rc.templates)
        tpl = mapper.validation(
            f = self.feature_info,
            linkpath = rc.templates.relpath().replace('\\','/'),
            overlap = overlap,
            p = self.performance,
            rc = rc,
            timestamp = statusfile.timestamp,
        ).respond(FileTransaction(rc.report_index, "w"))
