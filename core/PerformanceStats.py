"""Calculates performance statistics given the scores of the positive and
negative citations"""

from __future__ import division
import numpy as nx


                                     
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


class PerformanceStats:
    """Performance statistics calculation after cross validation.

    @group From constructor: pscores, nscores, alpha, P, N, A
    
    @ivar pscores: Increasing scores of positive articles.

    @ivar nscores: Increasing scores of negative articles.

    @ivar alpha: Balance of recall and precision in the F measure.

    @ivar P: Number of positive articles.

    @ivar N: Number of negative articles.

    @ivar A: Equal to L{P}+L{N}.
    
    
    @group From counts: uscores, vlen, PE, NE, TP, FN, FP, TP
    
    @ivar uscores: Unique scores in increasing order

    @ivar vlen: Length of performance vectors (= length of L{uscores})

    @ivar PE: Number of positives with each score in L{uscores}

    @ivar NE: Number of negatives with each score in L{uscores}
    
    @ivar TP, FN, FP, TN: Vectors for confusion matrix at each distinct threshold
    
    
    @group From make_ratio_vectors: TPR, FPR, PPV, FM, FMa
    
    @ivar TPR, FPR, PPV, FM, FMa: Vectors of performance ratios at each
    distinct threshold.
    
    
    @group From make_curve_areas: ROC_area,PR_area
    
    @ivar ROC_area: Area under ROC curve.
    
    @ivar PR_area: Aread under precision-recall curve.
    
    
    @ivar bep_index, breakeven: Breakeven point (where precision=recall), from
    L{find_breakeven}.
    
    @ivar threshold_index, threshold: Tuned threshold and its index, from
    L{maximise_fmeasure}.
    
    @ivar tuned: Tuned performance statistics, from L{get_tunedstats}.
    
    @ivar W, W_stderr: Better area under ROC curve, from L{roc_error}.
    
    @ivar AvPrec: Averaged precision, from L{averaged_precision}.
    """


    def __init__(self, pscores, nscores, alpha):
        """Constructor - parameters correspond to instance variables.

        @note: Sorted copies are made of L{pscores} and L{nscores}.
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
        s.make_confusion_matrix()
        s.make_ratio_vectors(s.alpha)
        s.make_curve_areas()
        s.roc_error()
        s.averaged_precision()
        s.maximise_fmeasure()
        s.find_breakeven()
        s.get_tunedstats()


    def make_confusion_matrix(self):
        """Calculates confusion matrix counts by iterating over pscores
        
        As a side effects, sets L{uscores}, L{vlen}, L{PE}, L{NE}
        
        @return: L{TP}, L{TN}, L{FP}, L{FN}
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
            ncount = TN # Start at TN and subtract it later
            while (ncount < s.N) and s.nscores[ncount] == threshold:
                ncount += 1
            s.NE[idx] = ncount-TN
            
            s.TP[idx] = TP
            s.TN[idx] = TN
            s.FP[idx] = FP
            s.FN[idx] = FN
        return s.TP, s.TN, s.FP, s.FN


    def make_ratio_vectors(self, alpha):
        """Calculate performance using vector algebra

        @param alpha: Weight of precision in calculating FMa
        
        @return: L{TPR}, L{FPR}, L{PPV}, L{FM}, L{FMa}
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


    def make_curve_areas(self):
        """Calculate areas under ROC and precision-recall curves
        
        Uses trapz(y, x). TPR is decreasing as threshold climbs, so vectors
        have to be reversed.
        
        This method underestimates ROC areas because boundary points (0,0) and
        (1,1) usually are not present in the data. Better to use L{roc_error}
        which does not have that problem.
        
        @return: L{ROC_area}, L{PR_area}"""
        from scipy.integrate import trapz
        s = self
        s.ROC_area = trapz(s.TPR[::-1], s.FPR[::-1])
        s.PR_area = trapz(s.PPV[::-1], s.TPR[::-1])
        return s.ROC_area, s.PR_area


    def mergescores(self):
        """Merged the contents of pscores and nscores in a single pass.
        
        Expects L{nscores} and L{pscores} in increasing order of score.        
        
        @return: Iterator over (score, relevance) in decreasing order of score.
        Relevance is True for members of pscores, and False for members of
        nscores. """
        s = self
        p_idx = s.P-1
        n_idx = s.N-1
        while p_idx >= 0 or n_idx >= 0:
            if p_idx >= 0 and \
            (n_idx < 0 or s.pscores[p_idx] >= s.nscores[n_idx]):
                yield s.pscores[p_idx], True
                p_idx -= 1
            elif n_idx >= 0 and \
            (p_idx < 0 or s.nscores[n_idx] > s.pscores[p_idx]):
                yield s.nscores[n_idx], False
                n_idx -= 1


    def averaged_precision(self):
        """Average the precision over each point of recall
        
        @return: L{AvPrec}, the precision averaged over each point where a
        relevant document is returned"""
        AvPrec = 0.0
        TP = 0
        FP = 0
        for score, relevant in self.mergescores():
            if relevant:
                TP += 1
                AvPrec += TP/(TP+FP)
            else:
                FP += 1
        self.AvPrec = AvPrec/TP
        return self.AvPrec


    def roc_error(self):
        """Area under ROC and its standard error
        
        Uses method of Hanley1982 to calculate standard error on the Wilcoxon
        statistic W, which corresponds to the area under the ROC by trapezoidal
        rule.
        
        @note: The vectors r1 .. r7 correspond to rows of Table II in
        Hanley1982.

        @return: L{W}, L{W_stderr}
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
        #print W, Q1, Q2, W_stderr
        s.W = W
        s.W_stderr = W_stderr
        return W, W_stderr


    def find_breakeven(self):
        """Calculate break-even, where precision equals recall.
        
        @return: L{bep_index}, L{breakeven} - index into pscores,
        and the recall/precision of the break-even point.
        """
        s = self
        diff = nx.absolute(nx.subtract(s.TPR, s.PPV))
        s.bep_index = nx.nonzero(diff == nx.min(diff))[0][0]
        s.breakeven = 0.5*(s.TPR[s.bep_index]+s.PPV[s.bep_index])
        return s.bep_index, s.breakeven


    def maximise_fmeasure(self):
        """Point of maximum F measure
        
        @return: L{threshold}, L{threshold_index}, the threshold and its index
        into uscores."""
        s = self
        max_FMa = nx.max(s.FMa)
        s.threshold_index = nx.nonzero(s.FMa == max_FMa)[0][0]
        s.threshold = self.uscores[s.threshold_index]
        return s.threshold, s.threshold_index


    def get_tunedstats(self):
        """Performance at the chosen threshold (usually the point of maximum F
        measure).
    
        @return: Storage object with these keys::
          P, N, A, T, F      (summary of input)
          TP, FP, TN, FN     (confusion matrix)
          TPR, FNR, TNR, FPR (ratios)
          PPV, NPV           (ratios)
          accuracy           (T/A)
          enrichment         (precision/prevalence)
          error              (F/A)
          fmeasure           (harmonic mean of TPR and PPV [alpha=0.5])
          fmeasure_alpha     (alpha-weighted F measure [alpha!=0.5])
          fmeasure_max       (maximum of standard F measure [alpha=0.5])
          precision          (PPV)
          prevalence         (P/A)
          recall             (TPR)
          specificity        (TNR)
          fp_tp_ratio        (FP/TP)
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
        TPR, FNR, TNR, FPR, PPV, NPV, FDR = 0, 0, 0, 0, 0, 0, 0
        if TP + FN != 0:
            TPR = TP / (TP + FN) # TPR=TP/P = sensitivity = recall
            FNR = FN / (TP + FN) # FNR=FN/P = 1-TP/P = 1-sensitivity = 1-recall
        if TN + FP != 0:
            TNR = TN / (TN + FP) # TNR=TN/N = specificity
            FPR = FP / (TN + FP) # FPR=FP/N = 1 - TN/N = 1-specificity
        if TP + FP != 0:
            PPV = TP / (TP + FP) # PPV=precision
            FDR = FP / (TP + FP) # FDR=1-precision
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
        # Return local variables in a Storage object
        from mscanner.core.Storage import Storage
        self.tuned = Storage(locals())
        del self.tuned.self
        return self.tuned