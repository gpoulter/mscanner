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
    """Cross-validated calculation of article scores
    
    Constructor arguments are also saved as instance variables.
    
    From calling validate():
      @ivar pscores: Scores of positive articles after validation
      @ivar nscores: Scores of negative articles after validation
    
    From calling getPerformance():
      @ivar performance: PerformanceStats object
    """

    def __init__(
        self, featmap, featdb, pos, neg, nfold,
        pseudocount = 0.01,
        alpha = 0.5,
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
        @param alpha: For maximising alpha-weighted F-Measure
        @param genedrug_articles: Set of doc ids with gene-drug co-occurrences
        @param mask: Features to exclude for scoring purposes
        @param norandom: If True forgoes randomisation in cross-validation for test purposes
        
        @note: If pseudocount is 0/False/None, we instead use the background
        medline frequency (n.b. not negative corpus frequency) of each 
        feature as its individual pseudocount.
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
        self.alpha = alpha
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
        """Carry out validation. 
        
        @return: self.pscores and self.nscores
        
        @note: self.nfolds==0 is taken to mean leave-out-one validation.
        """
        if self.nfold:
            return self.crossValidate()
        else:
            return self.leaveOutOneValidate()
        
    def getPerformance(self):
        """Return PerformanceStats object for this run 
        
        @note: Only call after validate()
        
        @return: self.performance
        """
        self.performance = PerformanceStats(self.pscores, self.nscores, self.alpha)
        return self.performance

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

        @return: self.pscores and self.nscores
        """
        pdocs = len(self.pos)
        ndocs = len(self.neg)
        log.info("%d pos and %d neg articles", pdocs, ndocs)
        if self.randomise:
            nx.random.shuffle(self.pos)
            nx.random.shuffle(self.neg)
        pstarts, psizes = self.partitionSizes(pdocs, self.nfold)
        nstarts, nsizes = self.partitionSizes(ndocs, self.nfold)
        self.pscores = nx.zeros(pdocs, nx.float32)
        self.nscores = nx.zeros(ndocs, nx.float32)
        pcounts = countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = countFeatures(self.numfeats, self.featdb, self.neg)
        for fold, (pstart,psize,nstart,nsize) in \
            enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            statusfile.update(fold, self.nfold)
            log.debug("pstart = %d, psize = %s; nstart = %d, nsize = %d", 
                      pstart, psize, nstart, nsize)
            # Modifiy the feature counts by subtracting out the test fold
            temp_pcounts = pcounts - countFeatures(
                self.numfeats, self.featdb, self.pos[pstart:pstart+psize])
            temp_ncounts = ncounts - countFeatures(
                self.numfeats, self.featdb, self.neg[nstart:nstart+nsize])
            # Calculate the resulting feature scores
            termscores, pfreqs, nfreqs = calculateFeatureScores(
                temp_pcounts, temp_ncounts, pdocs-psize, ndocs-nsize,
                self.pseudocount, self.mask, self.daniel)
            # Calculate the article scores for the test fold
            self.pscores[pstart:pstart+psize] = [
                nx.sum(termscores[self.featdb[d]]) for d in self.pos[pstart:pstart+psize]]
            self.nscores[nstart:nstart+nsize] = [
                nx.sum(termscores[self.featdb[d]]) for d in self.neg[nstart:nstart+nsize]]
        statusfile.update(None, self.nfold)
        return self.pscores, self.nscores
    
    def leaveOutOneValidate(self):
        """Carries out leave-out-one validation, returning the resulting scores.
        
        @note: Does not support 'daniel' scoring method
        
        @return: self.pscores and self.nscores
        """
        pcounts = countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = countFeatures(self.numfeats, self.featdb, self.neg)
        self.pscores = nx.zeros(len(self.pos), nx.float32)
        self.nscores = nx.zeros(len(self.neg), nx.float32)
        pdocs = len(self.pos)
        ndocs = len(self.neg)
        mask = self.mask
        if isinstance(self.pseudocount, float):
            ps = nx.zeros(self.numfeats, nx.float32) + self.pseudocount
        else:
            ps = self.pseudocount
        marker = 0
        for idx, doc in enumerate(self.pos):
            if idx == marker:
                statusfile.update(marker, pdocs+ndocs)
                marker += int((pdocs+ndocs)/20)
            f = [fid for fid in self.featdb[doc] if not mask or not mask[fid]]
            self.pscores[idx] = nx.sum(nx.log(((
                pcounts[f]-1+ps[f])/(pdocs-1+2*ps[f]))/((ncounts[f]+ps[f])/(ndocs+2*ps[f]))))
        for idx, doc in enumerate(self.neg):
            if pdocs+idx == marker:
                statfile(marker, pdocs+ndocs)
                marker += int((pdocs+ndocs)/20)
            f = [fid for fid in self.featdb[doc] if not mask or not mask[fid]]
            self.nscores[idx] = nx.sum(nx.log(((
                pcounts[f]+ps[f])/(pdocs+2*ps[f]))/((ncounts[f]-1+ps[f])/(ndocs-1+2*ps[f]))))
        statusfile.update(None, pdocs+ndocs)
        return self.pscores, self.nscores

class PerformanceStats:
    """Calculates and stores performance statistics.

    @ivar alpha, pscores, nscores: Input parameters
    
    @ivar P, N, A: Number of positive, negative articles (A=P+N)
    
    @ivar TP, FN, FP, FN, TPR, FPR, PPV, FM, FMa: Performance vectors
    
    @ivar ROC_area, PR_area: Integral statistics
    
    @ivar bep_index, breakeven: Index and value of breakeven point (where
    precision=recall)
    
    @ivar threshold_index, threshold: Index of and score value of the
    tuned threshold (cutoff for classification citations positive).
    
    @ivar tuned: Tuned performance statistics (see tunedStatistics)
    """

    def __init__(self, pscores, nscores, alpha):
        """Initialies the performance statistics. Parameters are kept
        as instance variables.

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @param alpha: Alpha parameter to balance recall and
        precision in FM_alpha.  Defaults to 0.5, producing harmonic
        mean of recall and precision (F-Measure)
        """
        s = self
        s.alpha = alpha
        s.pscores = pscores.copy()
        s.nscores = nscores.copy()
        s.pscores.sort()
        s.nscores.sort()
        s.P = len(s.pscores)
        s.N = len(s.nscores)
        s.A = s.P + s.N
        s.makeCountVectors()
        s.makeRatioVectors()
        s.makeCurveAreas()
        s.maxFMeasurePoint()
        s.breakEvenPoint()
        s.tunedStatistics()

    def makeCountVectors(self):
        """Calculates confusion matrix counts by iterating over pscores
        
        @return: self.(TP, TN, FP, FN)
        """
        s = self
        z = nx.zeros(s.P, nx.float32)
        s.TP = z.copy() # true positives
        s.TN = z.copy() # true negatives
        s.FP = z.copy() # false positives
        s.FN = z.copy() # false negatives
        TN = 1
        for xi in xrange(s.P):
            threshold = s.pscores[xi]
            while (s.nscores[TN-1] < threshold) and (TN < s.N):
                TN += 1
            while (s.nscores[TN-1] >= threshold) and (TN > 0):
                TN -= 1
            # xi+1 is how many positives we called negative
            FN = xi
            # TP+FN = P
            TP = s.P - FN
            # TN+FP = N
            FP = s.N - TN
            s.TP[xi] = TP
            s.TN[xi] = TN
            s.FP[xi] = FP
            s.FN[xi] = FN
        return s.TP, s.TN, s.FP, s.FN

    def makeRatioVectors(self):
        """Calculate performance using vector algebra
        
        @return: self.(TPR, FPR, PPV, FM, FMa)
        """
        s = self
        z = nx.zeros(s.P, nx.float32)
        # TPR is recall
        s.TPR = s.TP / s.P
        # FPR is 1-specificity
        s.FPR = s.FP / s.N
        # PPV is precision
        s.PPV = s.TP / (s.TP + s.FP) 
        s.PPV[s.TP+s.FP == 0] = 1.0
        # FM is F-Measure
        s.FM = 2 * s.TPR * s.PPV / (s.TPR + s.PPV) 
        # FMa is the alpha-weighted F-Measures
        s.FMa = 1 / (s.alpha / s.PPV + (1 - s.alpha) / s.TPR)
        return s.TPR, s.FPR, s.PPV, s.FM, s.FMa

    def makeCurveAreas(self):
        """Calculate areas under ROC and precision-recall curves
        
        @return: self.(ROC_area, PR_area)
        
        @note: trapz expects y, x. TPR is decreasing, so we reverse the vectors
        to present it in standard increasing form. ROC_area is calculated as 1
        - remaining area, since TPR covers 0.0 to 1.0 but FPR does not, meaning
        we must use TPR on the x-axis instead of FPR. """
        s = self
        s.ROC_area = 1.0 - trapz(s.FPR[::-1], s.TPR[::-1])
        s.PR_area = trapz(s.PPV[::-1], s.TPR[::-1])
        return s.ROC_area, s.PR_area
    
    def breakEvenPoint(self):
        """Calculate the threshold resulting in the break-even point
        where precision equals recall.  Returns index into pscores,
        and the recall/precision of the break-even point.
        
        @return: self.(bep_index, breakeven)
        """
        s = self
        diff = nx.absolute(nx.subtract(s.TPR, s.PPV))
        s.bep_index = nx.nonzero(diff == nx.min(diff))[0][0]
        s.breakeven = 0.5*(s.TPR[s.bep_index]+s.PPV[s.bep_index])
        return s.bep_index, s.breakeven

    def maxFMeasurePoint(self):
        """Calculate the threshold which results in the highest
        alpha-weighted F-Measure.  Returns the index into pscores for
        the threshold, and the threshold maximising F-Measure.
        
        @return: self.(threshold_index, threshold)
        """
        s = self
        s.threshold_index = nx.nonzero(s.FMa == nx.max(s.FMa))[0][0]
        s.threshold = self.pscores[s.threshold_index]
        return s.threshold_index, s.threshold

    def tunedStatistics(self):
        """Return Storage object with performance statistics for the
        tuned threshold.
    
        @return: Storage object with these keys:
          P, N, A, T, F    (summary of input)
          TP, FP, TN, FN   (confusion matrix)
          TPR, FNR, TNR, FPR
          PPV, NPV
          accuracy         (T/A)
          enrichment       (precision/prevalence)
          error            (F/A)
          fmeasure         (harmonic mean of TPR and PPV [alpha=0.5])
          fmeasure_alpha   (alpha-weighted F measure [alpha!=0.5])
          fmeasure_max     (maximum of standard F measure [alpha=0.5])
          precision        (PPV)
          prevalence       (P/A)
          recall           (TPR)
          specificity      (TNR)
        """
        TP = self.TP[self.threshold_index]
        TN = self.TN[self.threshold_index]
        FP = self.FP[self.threshold_index]
        FN = self.FN[self.threshold_index]
        P = self.P
        N = self.N
        A = self.A
        T = TP + TN
        F = FP + FN
        TPR, FNR, TNR, FPR, PPV, NPV = 0, 0, 0, 0, 0, 0
        if TP + FN != 0:
            TPR = TP / (TP + FN) # TPR=TP/P = sensitivity = recall
            FNR = FN / (TP + FN) # FNR=FN/P = 1-TP/P = 1-sensitivity = 1-recall
        if TN + FP != 0:
            TNR = TN / (TN + FP) # TNR=TN/N = specificity
            FPR = FP / (TN + FP) # FPR=FP/N = 1 - TN/N = 1-specificity
        if TP + FP != 0:
            PPV = TP / (TP + FP) # PPV=precision
        if TN + FN != 0:
            NPV = TN / (TN + FN) # NPV
        accuracy, prevalence = 0, 0
        if A != 0:
            accuracy = (TP + TN) / A   # acc  = T/A
            prevalence = (TP + FN) / A # prev = P/A
        error = 1 - accuracy
        recall = TPR
        specificity = TNR
        precision = PPV
        fmeasure, fmeasure_alpha, fmeasure_max = 0, 0, 0
        if recall > 0 and precision > 0:
            fmeasure = 2 * recall * precision / (recall + precision)
            fmeasure_alpha = 1.0 / ( (self.alpha / precision) + 
                                     ((1 - self.alpha) / recall))
            fmeasure_max = max(self.FM)
        enrichment = 0
        if prevalence > 0:
            enrichment = precision / prevalence
        # Return local variables as members of an object
        self.tuned = Storage(locals())
        del self.tuned.self
        return self.tuned 
