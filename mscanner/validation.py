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

from itertools import chain, izip
import logging as log
import numpy as nx
from scipy.integrate import trapz

from mscanner import statusfile
from mscanner.scoring import countFeatures
from mscanner.storage import Storage
from mscanner.utils import selfupdate

class Validator:
    """Cross-validated calculation of article scores
    
    Passed through constructor:
        @ivar featdb: Mapping from doc id to list of feature ids
        @ivar featinfo: FeatureInfo instance for stuff about features
        @ivar positives: Array of positive PMIDs for validation
        @ivar negatives: Array of negative PMIDs for validation
        @ivar nfolds: Number of validation folds (0 for leave-out-one)
        @ivar alpha: For Validator.getPerformance() PerformanceStatistics
        @ivar postfilter: Check score threshold AND membership of postfilter
    
    Updated by validate():
        @ivar pscores: Scores of positive articles after validation
        @ivar nscores: Scores of negative articles after validation
    
    Updated by getPerformance():
        @ivar performance: PerformanceStats object
    """

    def __init__(self, 
        featdb,
        featinfo,
        positives,
        negatives,
        nfolds,
        alpha = 0.5,
        postfilter = None,
        ):
        pscores = None
        nscores = None
        selfupdate()

    def articleIsPositive(self, docid, score, threshold):
        """Classifies an article as positive or negative based on score
        threshold
        
        @note: if self.positivestfilter is a non-empty set, we additionally
        check that the article is a member before classifying it positive. """
        if self.positivestfilter and docid in self.positivestfilter:
            return score >= threshold
        else:
            return False

    def validate(self):
        """Carry out validation. 
        
        @return: self.pscores and self.nscores
        
        @note: self.nfolds==0 is taken to mean leave-out-one validation.
        """
        if self.nfolds:
            return self.crossValidate()
        else:
            return self.leaveOutOneValidate()
        
    def getPerformance(self):
        """Return PerformanceStats object for this run 
        
        @note: Only call after validate()
        
        @return: self.performance
        """
        return PerformanceStats(self.pscores, self.nscores, self.alpha)

    @staticmethod
    def partitionSizes(nitems, nparts):
        """Returns start indeces and lengths of partitions for cross-validation"""
        base, rem = divmod(nitems, nparts)
        sizes = base * nx.ones(nparts, nx.int32)
        sizes[:rem] += 1
        starts = nx.zeros(nparts, nx.int32)
        starts[1:] = nx.cumsum(sizes[:-1])
        return starts, sizes

    def crossValidate(self, randomise=True):
        """Perform n-fold validation and return the raw performance measures
        
        @param randomise: Randomise validation splits (set False for debugging)
        
        @return: self.pscores and self.nscores
        """
        pdocs = len(self.positives)
        ndocs = len(self.negatives)
        log.info("%d pos and %d neg articles", pdocs, ndocs)
        if randomise:
            nx.random.shuffle(self.positives)
            nx.random.shuffle(self.negatives)
        pstarts, psizes = self.partitionSizes(pdocs, self.nfolds)
        nstarts, nsizes = self.partitionSizes(ndocs, self.nfolds)
        self.pscores = nx.zeros(pdocs, nx.float32)
        self.nscores = nx.zeros(ndocs, nx.float32)
        pcounts = countFeatures(len(self.featinfo), self.featdb, self.positives)
        ncounts = countFeatures(len(self.featinfo), self.featdb, self.negatives)
        for fold, (pstart,psize,nstart,nsize) in \
            enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            statusfile.update(fold, self.nfolds)
            log.debug("pstart = %d, psize = %s; nstart = %d, nsize = %d", 
                      pstart, psize, nstart, nsize)
            # Get new feature scores
            termscores = self.featinfo.updateFeatureScores(
                pos_counts = pcounts - countFeatures(
                    len(self.featinfo), self.featdb, 
                    self.positives[pstart:pstart+psize]), 
                neg_counts = ncounts - countFeatures(
                    len(self.featinfo), self.featdb, 
                    self.negatives[nstart:nstart+nsize]),
                pdocs = pdocs-psize, 
                ndocs = ndocs-nsize)
            # Calculate the article scores for the test fold
            self.pscores[pstart:pstart+psize] = [
                nx.sum(termscores[self.featdb[d]]) for d in 
                self.positives[pstart:pstart+psize]]
            self.nscores[nstart:nstart+nsize] = [
                nx.sum(termscores[self.featdb[d]]) for d in 
                self.negatives[nstart:nstart+nsize]]
        statusfile.update(None, self.nfolds)
        return self.pscores, self.nscores
    
    def leaveOutOneValidate(self):
        """Carries out leave-out-one validation, returning the resulting scores.
        
        @note: Is only able to perform Bayesian pseudocount feature score
        calculation.
        
        @return: self.pscores and self.nscores
        """
        pcounts = countFeatures(len(self.featinfo), self.featdb, self.positives)
        ncounts = countFeatures(len(self.featinfo), self.featdb, self.negatives)
        self.pscores = nx.zeros(len(self.positives), nx.float32)
        self.nscores = nx.zeros(len(self.negatives), nx.float32)
        pdocs = len(self.positives)
        ndocs = len(self.negatives)
        mask = self.featinfo.mask
        if isinstance(self.featinfo.pseudocount, float):
            ps = nx.zeros(len(self.featinfo), nx.float32) + self.featinfo.pseudocount
        else:
            ps = self.featinfo.pseudocount
        marker = 0
        for idx, doc in enumerate(self.positives):
            if idx == marker:
                statusfile.update(marker, pdocs+ndocs)
                marker += int((pdocs+ndocs)/20)
            f = [fid for fid in self.featdb[doc] if not mask or not mask[fid]]
            self.pscores[idx] = nx.sum(nx.log(
                ((pcounts[f]-1+ps[f])/(pdocs-1+2*ps[f]))/
                ((ncounts[f]+ps[f])/(ndocs+2*ps[f]))))
        for idx, doc in enumerate(self.negatives):
            if pdocs+idx == marker:
                statfile(marker, pdocs+ndocs)
                marker += int((pdocs+ndocs)/20)
            f = [fid for fid in self.featdb[doc] if not mask or not mask[fid]]
            self.nscores[idx] = nx.sum(nx.log(
                ((pcounts[f]+ps[f])/(pdocs+2*ps[f]))/
                ((ncounts[f]-1+ps[f])/(ndocs-1+2*ps[f]))))
        statusfile.update(None, pdocs+ndocs)
        return self.pscores, self.nscores
    
class PerformanceStats:
    """Stores performance statistics based on cross-validated article scores.
    
    @note: Effectively a table of (pscores, nscores, TP, FN, FP, TN, TPR,
    FPR, PPV, FM, FMa) with some scalar attributes.  Potential to
    store it in PyTables.

    @ivar alpha, pscores, nscores: Input parameters
    
    @ivar P, N, A: Number of positive, negative articles (A=P+N)
    
    @ivar TP, FN, FP, FN, TPR, FPR, PPV, FM, FMa: Performance vectors
    
    @ivar ROC_area, PR_area: Integral statistics
    
    @ivar bep_index, breakeven: Index and value of breakeven point (where
    precision=recall)
    
    @ivar threshold_index, threshold: Index of and score value of the
    tuned threshold (cutoff for classification citations positive).
    
    Set by tunedStatistics():
    
    @ivar tuned: Tuned performance statistics
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
          fp_tp_ratio      (FP/TP)
        """
        TP = int(self.TP[self.threshold_index])
        TN = int(self.TN[self.threshold_index])
        FP = int(self.FP[self.threshold_index])
        FN = int(self.FN[self.threshold_index])
        P = self.P  # P = TP + FN
        N = self.N  # N = TN + FP
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
        accuracy = T / A if A != 0 else 0
        prevalence = P / A if A != 0 else 0
        error = 1 - accuracy
        recall = TPR
        specificity = TNR
        precision = PPV
        fp_tp_ratio = FP/TP if TP != 0 else 0
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
