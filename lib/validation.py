"""Cross-validation and performance statistics module

@author: Graham Poulter
                                   

Validator -- Perform cross-validation and output performance statistics

"""

#Standard
from __future__ import division
import warnings
warnings.simplefilter('ignore', FutureWarning)
from itertools import chain, izip
import codecs
import cPickle
import logging as log
import math
import random
import time
import sys
#Local
from Gnuplot import Gnuplot
import numpy as nx
from scipy.integrate import trapz
from path import path
#MScanner
import article
import scoring
import templates
import plotting

class Validator:
    """Cross-validated calculation of article scores"""

    def __init__(
        self,
        numfeats,
        featdb,
        pos,
        neg,
        nfold,
        pseudocount = 0.01,
        daniel = False,
        genedrug_articles = None,
        mask = None,
        randomise = True,
        ):
        """Initialise validator
        @param numfeats: Number of distinct features (thus the length of the feature score vector)
        @param featdb: Mapping of doc id to list of feature ids
        @param pos: Array of positive doc ids
        @param neg: Array of negative doc ids
        @param nfold: Number of folds in cross-validation
        @param pseudocount, daniel: Parameters for scoring.calculateFeatureScores
        @param genedrug_articles: Set of doc ids with gene-drug co-occurrences
        @param mask: Features to exclude for scoring purposes
        @param norandom: If True forgoes randomisation in cross-validation for test purposes
        """
        self.numfeats = numfeats
        self.featdb = featdb
        self.pos = pos
        self.neg = neg
        self.nfold = nfold
        self.pseudocount = pseudocount
        self.daniel = daniel
        self.genedrug_articles = genedrug_articles
        self.mask = mask
        self.randomise = randomise

    def articleIsPositive(self, docid, score, threshold):
        """Classifies an article as positive or negative based on score threshold"""
        if self.genedrug_articles is None or docid in self.genedrug_articles:
            return score >= threshold
        else:
            return False

    def validate(self, statfile=None):
        """Carry out validation.  0 folds is taken to mean leave-out-one validation"""
        if self.nfold == 0:
            return self.leaveOutOneValidate(statfile)
        else:
            return self.crossValidate(statfile)

    @staticmethod
    def partitionSizes(nitems, nparts):
        """Returns start indeces and lengths of partitions for cross-validation"""
        base, rem = divmod(nitems, nparts)
        sizes = base * nx.ones(nparts, nx.int32)
        sizes[:rem] += 1
        starts = nx.zeros(nparts, nx.int32)
        starts[1:] = nx.cumsum(sizes[:-1])
        return starts, sizes

    @staticmethod
    def moveToFront(data, start, size):
        """Swaps a portion of data to the front of the array"""
        tmp = data[:size].copy()
        data[:size] = data[start:start+size]
        data[start:start+size] = tmp

    def crossValidate(self, statfile=None):
        """Perform n-fold validation and return the raw performance measures
        """
        log.info("%d pos and %d neg articles", len(self.pos), len(self.neg))
        if self.randomise:
            random.shuffle(self.pos)
            random.shuffle(self.neg)
        pos = self.pos.copy()
        neg = self.neg.copy()
        pstarts, psizes = self.partitionSizes(len(pos), self.nfold)
        nstarts, nsizes = self.partitionSizes(len(neg), self.nfold)
        pscores = nx.zeros(len(pos), nx.float32)
        nscores = nx.zeros(len(neg), nx.float32)
        pcounts = article.countFeatures(self.numfeats, self.featdb, pos)
        ncounts = article.countFeatures(self.numfeats, self.featdb, neg)
        for fold, (pstart,psize,nstart,nsize) in enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            log.debug("Carrying out fold number %d", fold)
            if statfile is not None:
                statfile.update(fold+1, self.nfold)
            # Move the test fold to the front 
            self.moveToFront(pos, pstart, psize)
            self.moveToFront(neg, nstart, nsize)
            log.debug("pstart = %d, psize = %s, nstart = %d, nsize = %d", pstart, psize, nstart, nsize)
            #log.debug("POS %s %s", repr(pos[:psize]), repr(pos[psize:]))
            #log.debug("NEG %s %s", repr(neg[:nsize]), repr(neg[nsize:]))
            # Modifiy the feature counts by subtracting out the test fold
            temp_pcounts = pcounts - article.countFeatures(self.numfeats, self.featdb, pos[:psize])
            temp_ncounts = ncounts - article.countFeatures(self.numfeats, self.featdb, neg[:nsize])
            #old_pcounts = article.countFeatures(self.numfeats, self.featdb, pos[psize:])
            #old_ncounts = article.countFeatures(self.numfeats, self.featdb, neg[nsize:])
            #log.debug("PDIFF %s", repr(old_pcounts - temp_pcounts))
            #log.debug("NDIFF %s", repr(old_ncounts - temp_ncounts))
            # Calculate the resulting feature scores
            termscores, pfreqs, nfreqs = scoring.calculateFeatureScores(
                temp_pcounts, temp_ncounts, len(pos)-psize, len(neg)-nsize,
                self.pseudocount, self.mask, self.daniel)
            # Calculate the article scores for the test fold
            pscores[pstart:pstart+psize] = [nx.sum(termscores[self.featdb[d]]) for d in pos[:psize]]
            nscores[nstart:nstart+nsize] = [nx.sum(termscores[self.featdb[d]]) for d in neg[:nsize]]
        return pscores, nscores
    
    def leaveOutOneValidate(self, statfile=None):
        """Carries out leave-out-one validation, returning the resulting scores.

        @note: Does not support 'daniel' scoring method
        """
        pcounts = article.countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = article.countFeatures(self.numfeats, self.featdb, self.neg)
        pscores = nx.zeros(len(self.pos), nx.float32)
        nscores = nx.zeros(len(self.neg), nx.float32)
        pdocs = len(self.pos)
        ndocs = len(self.neg)
        mask = None
        ignore_unknowns = True
        if ignore_unknowns:
            mask = (pcounts<=1) & (ncounts<=1)
        else:
            mask = nx.zeros(self.numfeats, nx.bool)
        if self.mask is not None:
            mask |= self.mask
        ps = self.pseudocount
        marker = 0
        if not statfile:
            marker = pdocs+ndocs+1
        for idx, doc in enumerate(self.pos):
            if idx == marker:
                statfile.update(marker, pdocs+ndocs)
                marker += (pdocs+ndocs)/20
            pscores[idx] = sum(
                math.log(((pcounts[f]-1+ps)/(pdocs-1+2*ps))/((ncounts[f]+ps)/(ndocs+2*ps))) for
                f in self.featdb[doc] if not mask[f])
        for idx, doc in enumerate(self.neg):
            if pdocs+idx == marker:
                statfile.update(marker, pdocs+ndocs)
                marker += (pdocs+ndocs)/20
            nscores[idx] = sum(
                math.log(((pcounts[f]+ps)/(pdocs+2*ps))/((ncounts[f]-1+ps)/(ndocs-1+2*ps))) for
                f in self.featdb[doc] if not mask[f])
        return pscores, nscores

def classDictionary(dictionary, exclude):
    """Converts a dictionary to a class with attributes of the dictionary keys,
    deleting those attributes listed in exclude."""
    class Empty: pass
    result = Empty()
    result.__dict__.update(dictionary)
    for attr in exclude:
        if hasattr(result, attr):
            delattr(result, attr)
    return result

class PerformanceStats:
    """Calculates and stores performance statistics.

    P, N, A, pscores, nscores, alpha,
    TP, FN, FP, FN, TPR, FPR, PPV, FM, FMa, 
    ROC_area, PR_area, 
    
    @ivar tuned: Object with tuned performance statistics (see tunedStatistics)
    """

    def __init__(self, pscores, nscores, alpha=0.5):
        """Initialies the performance statistics. Parameters are kept
        as instance variables.

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @param alpha: Alpha parameter to balance recall and
        precision in FM_alpha.  Defaults to 0.5, producing harmonic
        mean of recall and precision (F-Measure)
        """
        _ = self
        _.alpha = alpha
        _.pscores = pscores.copy()
        _.nscores = nscores.copy()
        _.pscores.sort()
        _.nscores.sort()
        _.P = len(_.pscores)
        _.N = len(_.nscores)
        _.A = _.P + _.N
        _.makeCountVectors()
        _.makeRatioVectors()
        _.makeCurveAreas()
        _.maxFMeasurePoint()
        _.breakEvenPoint()
        _.tunedStatistics()

    def makeCountVectors(self):
        """Calculates arrays of TP, TN, FP, FN by iterating over pscores"""
        _ = self
        z = nx.zeros(_.P, nx.float32)
        _.TP = z.copy() # true positives
        _.TN = z.copy() # true negatives
        _.FP = z.copy() # false positives
        _.FN = z.copy() # false negatives
        TN = 1
        for xi in xrange(_.P):
            threshold = _.pscores[xi]
            while (_.nscores[TN-1] < threshold) and (TN < _.N):
                TN += 1
            while (_.nscores[TN-1] >= threshold) and (TN > 0):
                TN -= 1
            # xi+1 is how many positives we called negative
            FN = xi
            # TP+FN = P
            TP = _.P - FN
            # TN+FP = N
            FP = _.N - TN
            _.TP[xi] = TP
            _.TN[xi] = TN
            _.FP[xi] = FP
            _.FN[xi] = FN
        return _.TP, _.TN, _.FP, _.FN

    def makeRatioVectors(self):
        """Calculate arrays of TPR, FPR, PPV, FM, FMa"""
        _ = self
        z = nx.zeros(_.P, nx.float32)
        # TPR is recall
        _.TPR = _.TP / _.P
        # FPR is 1-specificity
        _.FPR = _.FP / _.N
        # PPV is precision
        _.PPV = _.TP / (_.TP + _.FP) 
        _.PPV[_.TP+_.FP == 0] = 1.0
        # FM is F-Measure
        _.FM = 2 * _.TPR * _.PPV / (_.TPR + _.PPV) 
        # FMa is the alpha-weighted F-Measures
        _.FMa = 1 / (_.alpha / _.PPV + (1 - _.alpha) / _.TPR)
        return _.TPR, _.FPR, _.PPV, _.FM, _.FMa

    def makeCurveAreas(self):
        """Calculate areas under ROC and precision-recall curves"""
        _ = self
        # Note that trapz expects y, x TPR is decreasing, so we reverse the
        # vectors to present it in standard increasing form. ROC_area is
        # calculated as 1 - remaining area, since TPR covers 0.0 to 1.0 but FPR
        # does not, meaning we must use TPR on the x-axis instead of FPR.
        _.ROC_area = 1.0 - trapz(_.FPR[::-1], _.TPR[::-1])
        _.PR_area = trapz(_.PPV[::-1], _.TPR[::-1])
        return _.ROC_area, _.PR_area
    
    def breakEvenPoint(self):
        """Calculate the threshold resulting in the break-even point
        where precision equals recall.  Returns index into pscores,
        and the recall/precision of the break-even point.
        """
        _ = self
        diff = nx.absolute(nx.subtract(_.TPR, _.PPV))
        _.bep_index = nx.nonzero(diff == nx.min(diff))[0][0]
        _.breakeven = 0.5*(_.TPR[_.bep_index]+_.PPV[_.bep_index])
        return _.bep_index, _.breakeven

    def maxFMeasurePoint(self):
        """Calculate the threshold which results in the highest
        alpha-weighted F-Measure.  Returns the index into pscores for
        the threshold, and the threshold maximising F-Measure.
        """
        _ = self
        _.fmax_index = nx.nonzero(_.FMa == nx.max(_.FMa))[0][0]
        _.threshold = self.pscores[_.fmax_index]
        return _.fmax_index, _.threshold

    def tunedStatistics(self):
        """Object contains performance statistics for a particular tuned threshold.
    
        Instance variables include:
        threshold,
        TP, FP, TN, FN,
        P, N, A, T, F,
        TPR, FNR, TNR, FPR, PPV, NPV,
        accuracy, prevalence,
        recall, precision,
        fmeasure, fmeasure_alpha, fmeasure_max,
        enrichment
        """
        threshold = self.threshold
        TP = self.TP[self.fmax_index]
        TN = self.TN[self.fmax_index]
        FP = self.FP[self.fmax_index]
        FN = self.FN[self.fmax_index]
        P = self.P
        N = self.N
        A = self.A
        T = TP + TN
        F = FP + FN
        TPR, FNR, TNR, FPR, PPV, NPV = 0, 0, 0, 0, 0, 0
        if TP + FN != 0:
            TPR = TP / (TP + FN) # TPR = TP/P = sensitivity = recall
            FNR = FN / (TP + FN) # FNR = FN/P = 1 - TP/P = 1-sensitivity = 1-recall
        if TN + FP != 0:
            TNR = TN / (TN + FP) # TNR = TN/N = specificity
            FPR = FP / (TN + FP) # FPR = FP/N = 1 - TN/N = 1-specificity
        if TP + FP != 0:
            PPV = TP / (TP + FP) # PPV = precision
        if TN + FN != 0:
            NPV = TN / (TN + FN) # NPV
        accuracy, prevalence = 0, 0
        if A != 0:
            accuracy = (TP + TN) / A   # acc  = T/A
            prevalence = (TP + FN) / A # prev = P/A
        recall = TPR
        precision = PPV
        fmeasure, fmeasure_alpha, fmeasure_max = 0, 0, 0
        if recall > 0 and precision > 0:
            fmeasure = 2 * recall * precision / (recall + precision)
            fmeasure_alpha = 1.0 / (self.alpha / precision + (1 - self.alpha) / recall)
            fmeasure_max = max(self.FM)
        enrichment = 0
        if prevalence > 0:
            enrichment = precision / prevalence
        # Return local variables as members of an object
        self.tuned = classDictionary(locals(), ["self"])
        return self.tuned 

def report(pos, neg, pscores, nscores, featmap, featdb, configuration):
    """Write a full validation report
    
    @param pos: IDs of positive articles
    @param neg: IDs of negative articles
    
    @param pscores: Scores of positive articles
    @param nscores: Scores of negative articles

    @param featmap: Mapping feature ID <-> feature string
    @param featdb: Mapping doc ID -> list of feature IDs
    @param configuration: Configuration module to source parameters
    """
    c = configuration
    rd = c.reportdir
    p = PerformanceStats(pscores, nscores, c.alpha)
    gp = Gnuplot(debug=1)
    feature_info = scoring.FeatureScoreInfo(
        pos_counts = article.countFeatures(len(featmap), featdb, pos),
        neg_counts = article.countFeatures(len(featmap), featdb, neg),
        pdocs = len(pos),
        ndocs = len(neg),
        pseudocount = c.pseudocount,
        featmap = featmap,
        exclude_types = c.exclude_types,
        daniel = c.dodaniel
        )
    feature_info.writeFeatureScoresCSV(codecs.open(rd/c.term_scores_name,"wb","utf-8"))
    overlap, iX, iY = plotting.plotPDF(
        gp, rd/c.hist_img, p.pscores, p.nscores, p.tuned.threshold)
    plotting.plotROC(
        gp, rd/c.roc_img, p.FPR, p.TPR, p.tuned.FPR)
    plotting.plotPrecisionRecall(
        gp, rd/c.p_vs_r_img, p.TPR, p.PPV, p.tuned.TPR)
    plotting.plotPrecisionRecallFmeasure(
        gp, rd/c.pr_vs_score_img, p.pscores, p.TPR, p.PPV, p.FM, p.FMa, p.tuned.threshold)
    # Write main index file for validation output
    templates.validation.run(dict(
        time = time.strftime("%Y-%m-%d %H:%M:%S"),
        overlap = overlap,
        p = p,
        t = p.tuned,
        f = feature_info,
        c = configuration,
        ), outputfile=file(rd/c.index_file, "w"))
    c.stylesheet.copy(rd/"style.css")
