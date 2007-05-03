"""Cross-validation and performance statistics module

Validator -- Perform cross-validation for article scores
PerformanceStats -- Calculate statistics from article scores

                                   
"""

from __future__ import division

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

import warnings
warnings.simplefilter('ignore', FutureWarning)

import logging as log
import numpy as nx
from scipy.integrate import trapz

from mscanner import statusfile
from mscanner.scoring import calculateFeatureScores
from mscanner.storage import Storage
from mscanner.utils import countFeatures

class Validator:
    """Cross-validated calculation of article scores"""

    def __init__(
        self, featmap, featdb, pos, neg, nfold,
        pseudocount = 0.01,
        daniel = False,
        genedrug_articles = None,
        mask = None,
        randomise = True,
        ):
        """Initialise validator
        @param featmap: FeatureMapping instance (for fixed pseudocount only len()is used)
        @param featdb: Mapping of doc id to list of feature ids
        @param pos: Array of positive doc ids
        @param neg: Array of negative doc ids
        @param nfold: Number of folds in cross-validation
        @param pseudocount, daniel: Parameters for scoring.calculateFeatureScores
        @param genedrug_articles: Set of doc ids with gene-drug co-occurrences
        @param mask: Features to exclude for scoring purposes
        @param norandom: If True forgoes randomisation in cross-validation for test purposes
        """
        self.featmap = featmap
        self.numfeats = len(featmap)
        self.featdb = featdb
        self.pos = pos
        self.neg = neg
        self.nfold = nfold
        self.pseudocount = pseudocount
        if not pseudocount:
            self.pseudocount = nx.array(featmap.counts, nx.float32) / featmap.numdocs
        self.daniel = daniel
        self.genedrug_articles = genedrug_articles
        self.mask = mask
        self.randomise = randomise

    def articleIsPositive(self, docid, score, threshold):
        """Classifies an article as positive or negative based on score threshold"""
        if self.genedrug_articles and docid in self.genedrug_articles:
            return score >= threshold
        else:
            return False

    def validate(self):
        """Carry out validation.  0 folds is taken to mean leave-out-one validation"""
        if self.nfold:
            return self.crossValidate()
        else:
            return self.leaveOutOneValidate()

    @staticmethod
    def partitionSizes(nitems, nparts):
        """Returns start indeces and lengths of partitions for cross-validation"""
        base, rem = divmod(nitems, nparts)
        sizes = base * nx.ones(nparts, nx.int32)
        sizes[:rem] += 1
        starts = nx.zeros(nparts, nx.int32)
        starts[1:] = nx.cumsum(sizes[:-1])
        return starts, sizes

    def crossValidate(self):
        """Perform n-fold validation and return the raw performance measures
        """
        pdocs = len(self.pos)
        ndocs = len(self.neg)
        log.info("%d pos and %d neg articles", pdocs, ndocs)
        if self.randomise:
            nx.random.shuffle(self.pos)
            nx.random.shuffle(self.neg)
        pstarts, psizes = self.partitionSizes(pdocs, self.nfold)
        nstarts, nsizes = self.partitionSizes(ndocs, self.nfold)
        pscores = nx.zeros(pdocs, nx.float32)
        nscores = nx.zeros(ndocs, nx.float32)
        pcounts = countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = countFeatures(self.numfeats, self.featdb, self.neg)
        for fold, (pstart,psize,nstart,nsize) in enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            statusfile.update(fold, self.nfold)
            log.debug("pstart = %d, psize = %s; nstart = %d, nsize = %d", pstart, psize, nstart, nsize)
            # Modifiy the feature counts by subtracting out the test fold
            temp_pcounts = pcounts - countFeatures(self.numfeats, self.featdb, self.pos[pstart:pstart+psize])
            temp_ncounts = ncounts - countFeatures(self.numfeats, self.featdb, self.neg[nstart:nstart+nsize])
            # Calculate the resulting feature scores
            termscores, pfreqs, nfreqs = calculateFeatureScores(
                temp_pcounts, temp_ncounts, pdocs-psize, ndocs-nsize,
                self.pseudocount, self.mask, self.daniel)
            # Calculate the article scores for the test fold
            pscores[pstart:pstart+psize] = [nx.sum(termscores[self.featdb[d]]) for d in self.pos[pstart:pstart+psize]]
            nscores[nstart:nstart+nsize] = [nx.sum(termscores[self.featdb[d]]) for d in self.neg[nstart:nstart+nsize]]
        statusfile.update(None, self.nfold)
        return pscores, nscores
    
    def leaveOutOneValidate(self):
        """Carries out leave-out-one validation, returning the resulting scores.

        @note: Does not support 'daniel' scoring method
        """
        pcounts = countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = countFeatures(self.numfeats, self.featdb, self.neg)
        pscores = nx.zeros(len(self.pos), nx.float32)
        nscores = nx.zeros(len(self.neg), nx.float32)
        pdocs = len(self.pos)
        ndocs = len(self.neg)
        mask = None
        #mask = (pcounts<=1) & (ncounts<=1)
        if self.mask is not None:
            if mask is None:
                mask = self.mask
            else:
                mask |= self.mask
        if isinstance(self.pseudocount, float):
            ps = nx.zeros(self.numfeats, nx.float32) + self.pseudocount
        else:
            ps = self.pseudocount
        marker = 0
        for idx, doc in enumerate(self.pos):
            if idx == marker:
                statusfile.update(marker, pdocs+ndocs)
                marker += int((pdocs+ndocs)/20)
            f = [fid for fid in self.featdb[doc] if mask is None or not mask[fid]]
            pscores[idx] = nx.sum(nx.log(((pcounts[f]-1+ps[f])/(pdocs-1+2*ps[f]))/((ncounts[f]+ps[f])/(ndocs+2*ps[f]))))
        for idx, doc in enumerate(self.neg):
            if pdocs+idx == marker:
                statfile(marker, pdocs+ndocs)
                marker += int((pdocs+ndocs)/20)
            f = [fid for fid in self.featdb[doc] if mask is None or not mask[fid]]
            nscores[idx] = nx.sum(nx.log(((pcounts[f]+ps[f])/(pdocs+2*ps[f]))/((ncounts[f]-1+ps[f])/(ndocs-1+2*ps[f]))))
        statusfile.update(None, pdocs+ndocs)
        return pscores, nscores

class PerformanceStats:
    """Calculates and stores performance statistics.

    @ivar P, N, A, alpha, pscores, nscores: Input variables
    
    @ivar TP, FN, FP, FN, TPR, FPR, PPV, FM, FMa: Performance graph data
    
    @ivar ROC_area, PR_area: Integral statistics
    
    @ivar tuned: Tuned performance statistics (see tunedStatistics)
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
        self.tuned = Storage(locals())
        del self.tuned.self
        return self.tuned 
