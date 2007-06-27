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

from mscanner.scoring import countFeatures
from mscanner.support.utils import selfupdate

class Validator:
    """Cross-validated calculation of article scores
    
    Passed through constructor:
        @ivar featdb: Mapping from doc id to list of feature ids
        @ivar featinfo: FeatureInfo instance for stuff about features
        @ivar positives: Array of positive PMIDs for validation
        @ivar negatives: Array of negative PMIDs for validation
        @ivar nfolds: Number of validation folds (0 for leave-out-one)
        @ivar alpha: For Validator.getPerformance() PerformanceStatistics
        @ivar postfilter: Check membership of this as well as score threshold
    
    Updated by validate():
        @ivar pscores: Scores of positive articles after validation
        @ivar nscores: Scores of negative articles after validation
    """

    def __init__(self, 
        featdb,
        featinfo,
        positives,
        negatives,
        nfolds = 10,
        alpha = 0.5,
        postfilter = None,
        ):
        pscores = None
        nscores = None
        selfupdate()

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
        s = self
        pdocs = len(s.positives)
        ndocs = len(s.negatives)
        log.info("%d pos and %d neg articles", pdocs, ndocs)
        if randomise:
            nx.random.shuffle(s.positives)
            nx.random.shuffle(s.negatives)
        s.pstarts, s.psizes = s.partitionSizes(pdocs, s.nfolds)
        s.nstarts, s.nsizes = s.partitionSizes(ndocs, s.nfolds)
        s.pscores = nx.zeros(pdocs, nx.float32)
        s.nscores = nx.zeros(ndocs, nx.float32)
        pcounts = countFeatures(len(s.featinfo), s.featdb, s.positives)
        ncounts = countFeatures(len(s.featinfo), s.featdb, s.negatives)
        for fold, (pstart,psize,nstart,nsize) in \
            enumerate(zip(s.pstarts,s.psizes,s.nstarts,s.nsizes)):
            log.debug("Fold %d: pstart = %d, psize = %s; nstart = %d, nsize = %d", 
                      fold, pstart, psize, nstart, nsize)
            # Get new feature scores
            termscores = s.featinfo.updateFeatureScores(
                pos_counts = pcounts - countFeatures(
                    len(s.featinfo), s.featdb, 
                    s.positives[pstart:pstart+psize]), 
                neg_counts = ncounts - countFeatures(
                    len(s.featinfo), s.featdb, 
                    s.negatives[nstart:nstart+nsize]),
                pdocs = pdocs-psize, 
                ndocs = ndocs-nsize)
            # Calculate the article scores for the test fold
            s.pscores[pstart:pstart+psize] = [
                nx.sum(termscores[s.featdb[d]]) for d in 
                s.positives[pstart:pstart+psize]]
            s.nscores[nstart:nstart+nsize] = [
                nx.sum(termscores[s.featdb[d]]) for d in 
                s.negatives[nstart:nstart+nsize]]
        return s.pscores, s.nscores
    
    def leaveOutOneValidate(self):
        """Carries out leave-out-one validation, returning the resulting scores.
        
        @note: Updates self.pscores, self.nscores
        
        @note: Feature scores by Bayesian pseudocount only - no other methods.
        
        @return: self.pscores and self.nscores
        """
        # Set up base feature scores
        pcounts = countFeatures(len(self.featinfo), self.featdb, self.positives)
        ncounts = countFeatures(len(self.featinfo), self.featdb, self.negatives)
        self.pscores = nx.zeros(len(self.positives), nx.float32)
        self.nscores = nx.zeros(len(self.negatives), nx.float32)
        pdocs = len(self.positives)
        ndocs = len(self.negatives)
        mask = self.featinfo.mask
        # Set up pseudocount
        if isinstance(self.featinfo.pseudocount, float):
            ps = nx.zeros(len(self.featinfo), nx.float32) + self.featinfo.pseudocount
        else:
            ps = self.featinfo.pseudocount
        marker = 0
        # Discount this article in feature score calculations
        def score_article(pmid, p_mod, n_mod):
            f = [fid for fid in self.featdb[doc] if not mask or not mask[fid]]
            return nx.sum(nx.log(
                ((pcounts[f]+p_mod+ps[f])/(pdocs+p_mod+2*ps[f]))/
                ((ncounts[f]+n_mod+ps[f])/(ndocs+n_mod+2*ps[f]))))
        # Get scores for positive articles
        for idx, doc in enumerate(self.positives):
            self.pscores[idx] = score_article(doc, -1, 0)
        # Get scores for negative articles
        for idx, doc in enumerate(self.negatives):
            self.nscores[idx] = score_article(doc, 0, -1)
        return self.pscores, self.nscores
    
class PerformanceStats:
    """Stores performance statistics based on cross-validated article scores.
    
    Set by constructor:
    
        @ivar pscores: Scores of positive articles in increasing order

        @ivar nscores: Scores of negative articles in increasing order
    
        @ivar alpha: Balance of recall and precision in FM_alpha. Defaults to 0.5,
        producing harmonic mean of recall and precision (F-Measure)
    
        @ivar P, N, A: Number of positive, negative articles (A=P+N)
    
    Set by getCountVectors():
    
        @ivar uscores: Unique scores in increasing order
    
        @ivar vlen: Length of performance vectors (= length of uscores)

        @ivar PE, NE: Number of positives, negatives having each score in uscores.
    
        @ivar TP, FN, FP, FN: Confusion matrix for threshold at each score in
        uscores
        
    Set by getRatioVectors():
    
        @ivar TPR, FPR, PPV, FM, FMa: Performance vectors,
        with a value at each distinct score level.
        
    Set by getCurveAreas():
    
        @ivar ROC_area: Area under ROC curve
        
        @ivar PR_area: Aread under precision-recall curve
        
    Set by getBreakEvenPoint():
    
        @ivar bep_index, breakeven: Index and value of breakeven point (where
        precision=recall)
        
    Set by getMaxFMeasurePoint():
    
        @ivar threshold_index, threshold: Index of and score value of the
        tuned threshold (cutoff for classification citations positive).
    
    Set by tunedStatistics():
    
        @ivar tuned: Tuned performance statistics
    """

    def __init__(self, pscores, nscores, alpha):
        """Initialies the performance statistics. Parameters are 
        defined under instance variables."""
        s = self
        s.alpha = alpha
        s.pscores = pscores.copy()
        s.nscores = nscores.copy()
        s.pscores.sort()
        s.nscores.sort()
        s.P = len(s.pscores)
        s.N = len(s.nscores)
        s.A = s.P + s.N
        s.getCountVectors()
        s.getRatioVectors(s.alpha)
        s.getCurveAreas()
        s.getROCError()
        s.getMaxFMeasurePoint()
        s.getBreakEvenPoint()
        s.getTunedStatistics()
        
    def getCountVectors(self):
        """Calculates confusion matrix counts by iterating over pscores
        
        @return: self.(TP, TN, FP, FN)
        """
        s = self
        s.uscores = nx.unique(nx.concatenate((s.pscores,s.nscores)))
        s.vlen = len(s.uscores)
        s.PE = nx.zeros(s.vlen, nx.float32) # positives with given score
        s.NE = nx.zeros(s.vlen, nx.float32) # negatives with given score
        s.TP = nx.zeros(s.vlen, nx.float32) # true positives
        s.TN = nx.zeros(s.vlen, nx.float32) # true negatives
        s.FP = nx.zeros(s.vlen, nx.float32) # false positives
        s.FN = nx.zeros(s.vlen, nx.float32) # false negatives
        TN = 0
        FN = 0
        for idx, threshold in enumerate(s.uscores):

            # Classify positives scoring < threshold as negative
            # (look up score for the next article to classify)
            while (FN < s.P) and (s.pscores[FN] < threshold):
                FN += 1
            TP = s.P - FN # TP+FN=P

            # pcount-FN = number of positives having threshold score
            pcount = FN # Start at FN, subtract it later
            while (pcount < s.P) and s.pscores[pcount] == threshold:
                pcount += 1
            s.PE[idx] = pcount-FN
            
            # Classify negatives scoring < threshold as negative
            while (TN < s.N) and (s.nscores[TN] < threshold):
                TN += 1
            FP = s.N - TN  # TN+FP=N

            # ncount-TN = number of negatives having threshold score
            ncount = TN 
            while (ncount < s.N) and s.nscores[ncount] == threshold:
                ncount += 1
            s.NE[idx] = ncount-TN
            
            s.TP[idx] = TP
            s.TN[idx] = TN
            s.FP[idx] = FP
            s.FN[idx] = FN
        return s.TP, s.TN, s.FP, s.FN

    def getRatioVectors(self, alpha):
        """Calculate performance using vector algebra

        @param alpha: Weight of precision in calculating FMa
        
        @return: self.(TPR, FPR, PPV, FM, FMa)
        """
        s = self
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
        s.FMa = 1 / (alpha / s.PPV + (1 - alpha) / s.TPR)
        return s.TPR, s.FPR, s.PPV, s.FM, s.FMa

    def getCurveAreas(self):
        """Calculate areas under ROC and precision-recall curves
        
        @return: self.ROC_area, self.PR_area
        
        @note: Uses trapz(y, x). TPR is decreasing as threshold climbs, so
        vectors have to be reversed.
        
        @note: Underestimates ROC areas because boundary points (0,0) and (1,1)
        usually are not present in the data. Better to use W (from
        getROCError) which does not have that problem. """
        from scipy.integrate import trapz
        s = self
        s.ROC_area = trapz(s.TPR[::-1], s.FPR[::-1])
        s.PR_area = trapz(s.PPV[::-1], s.TPR[::-1])
        return s.ROC_area, s.PR_area
    
    def getROCError(self):
        """Use Hanley1982 to calculate standard error on the Wilcoxon
        statistic W, which corresponds to the area under the ROC
        by trapezoidal rule. 
        
        @note: The vectors r1 .. r7 correspond to rows of Table II in
        Hanley1982.
        
        @note: W, the Wilcoxon statistic, is equal to the trapezoidal
        rule area under the ROC.
        """
        s = self
        # r1 is number of negatives with each score,
        # r2 is number of positives rated higher than each score
        # r3 is number of positives with each score
        # r4 is number of negatives rated lower than each score
        r1 = s.NE
        r2 = s.TP - s.PE
        r3 = s.PE
        r4 = s.TN
        r5 = r1 * r2 + 0.5 * r1 * r3
        r6 = r3 * (r4**2 + r4*r1 + (r1**2)/3)
        r7 = r1 * (r2**2 + r2*r3 + (r3**2)/3)
        N = float(s.N)
        P = float(s.P)
        W = r5.sum() / (N*P)
        Q2 = r6.sum() / (P * N**2)
        Q1 = r7.sum() / (N * P**2)
        W_stderr = nx.sqrt((W*(1-W)+(P-1)*(Q1-W**2)+(N-1)*(Q2-W**2))/(P*N))
        print W, Q1, Q2, W_stderr
        s.W = W
        s.W_stderr = W_stderr
        return W, W_stderr
    
    def getBreakEvenPoint(self):
        """Calculate the threshold resulting in the break-even point
        where precision equals recall.  Returns index into pscores,
        and the recall/precision of the break-even point.
        
        @return: bep_index, breakeven
        """
        s = self
        diff = nx.absolute(nx.subtract(s.TPR, s.PPV))
        s.bep_index = nx.nonzero(diff == nx.min(diff))[0][0]
        s.breakeven = 0.5*(s.TPR[s.bep_index]+s.PPV[s.bep_index])
        return s.bep_index, s.breakeven

    def getMaxFMeasurePoint(self):
        """Return the threshold and its index into uscores that gives the
        greates value for alpha-weighted F-Measure."""
        s = self
        max_FMa = nx.max(s.FMa)
        s.threshold_index = nx.nonzero(s.FMa == max_FMa)[0][0]
        s.threshold = self.uscores[s.threshold_index]
        return s.threshold, s.threshold_index

    def getTunedStatistics(self):
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
        from mscanner.support.storage import Storage
        self.tuned = Storage(locals())
        del self.tuned.self
        return self.tuned
