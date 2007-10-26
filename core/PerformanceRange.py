"""Given a threshold, find the minimum and maximum for the precision, recall, etc.
across the validation folds."""

from __future__ import division
import numpy as nx
from mscanner.core.Validator import CrossValidator
from mscanner import update


                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"



class PerformanceRange:
    """Given a threshold, find the minimum and maximum for the precision,
    recall across the validation folds.

    @group pscores, nscores, nfolds, threshold: Passed via constructor
    
    @ivar pscores: Unsorted scores of positive documents.
    @ivar nscores: Unsorted scores of negative documents.
    @ivar nfolds: Number of cross validation folds
    @ivar threshold: Documents scoring above this are predicted positive.
    
    @ivar TP: Vector with TP at threshold for each fold
    @ivar FP: Vector with FP at threshold for each fold
    @ivar TN: Vector with TN at threshold for each fold
    @ivar FN: Vector with FN at threshold for each fold

    @ivar precision: (min, max) of precision across all folds
    @ivar recall: (min, max) of recall across all folds
    @ivar fmeasure: (min, max) of F measure across all folds
    """
    
    def __init__(self, pscores, nscores, nfolds, threshold):
        """Parameters correspond to instance variables"""
        update(self, locals())
        self.do_confusion_vectors()
        self.do_performance_range()


    def do_confusion_vectors(self):
        """Finds TP, TN, FP, FN at the threshold over each validation fold"""
        for vname in ["TP", "TN", "FP", "FN"]:
            setattr(self, vname, nx.zeros(self.nfolds, nx.float32))
        pstarts, psizes = CrossValidator.make_partitions(
            len(self.pscores), self.nfolds)
        nstarts, nsizes = CrossValidator.make_partitions(
            len(self.nscores), self.nfolds)
        for fold, (pstart,psize,nstart,nsize) in \
            enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            self.do_confusion_single(
                fold, 
                self.pscores[pstart:pstart+psize],
                self.nscores[nstart:nstart+nsize])


    def do_confusion_single(self, fold, pos, neg):
        """Find TP, TN, FP, FN at threshold inside a single validation fold"""
        # Find False Negatives and True Positives
        pos = nx.array(pos)
        pos.sort()
        P = len(pos)
        FN = 0
        while (FN < P) and (pos[FN] < self.threshold):
            FN += 1
        self.FN[fold] = FN
        self.TP[fold] = P - FN # TP+FN=P
        # Find True Negatives and False Positives
        neg = nx.array(neg)
        neg.sort()
        N = len(neg)
        TN = 0
        while (TN < N) and (neg[TN] < self.threshold):
            TN += 1
        self.TN[fold] = TN
        self.FP[fold] = N - TN # TN+FP=N


    def do_performance_range(self):
        """Finds (min,max) of precision, etc., using the TP/TN/FP/FN vectors
        over the folds."""
        for vname in ["precision", "recall", "fmeasure"]:
            setattr(self, vname, (1.0, 0.0))
        for TP, FP, TN, FN in zip(self.TP, self.FP, self.TN, self.FN):
            prec = (TP/(TP+FP)) if (TP+FP>0) else 0
            rec = (TP/(TP+FN))
            F = 2 * prec * rec / (prec + rec) if (prec+rec>0) else 0
            self._minimax("precision", prec)
            self._minimax("recall", rec)
            self._minimax("fmeasure", F)
    

    def _minimax(self, varname, value):
        """Given the name of a (min,max) attribute, update it
        with a new value if that value extends the range.
        
        @param varname: Name of attribute to update
        @param value: If value is outside the range
        """
        cmin, cmax = getattr(self, varname)
        if value < cmin:
            cmin = value
        if value > cmax:
            cmax = value
        setattr(self, varname, (cmin, cmax))