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
import numpy
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
        pseudocount = 0.1,
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
            return self.leave_out_one_validate(statfile)
        else:
            #return self.cross_validate(statfile)
            return self.cross_validate_new(statfile)

    @staticmethod
    def partitionSizes(nitems, nparts):
        """Returns start indeces and lengths of partitions for cross-validation"""
        base, rem = divmod(nitems, nparts)
        sizes = base * numpy.ones(nparts, dtype=numpy.int32)
        sizes[:rem] += 1
        starts = numpy.zeros(nparts, dtype=numpy.int32)
        starts[1:] = numpy.cumsum(sizes[:-1])
        return starts, sizes

    def cross_validate_new(self, statfile=None):
        """Perform n-fold validation and return the raw performance measures

        Calculates term scores a la leave_one_out_validate
        """
        def moveToFront(data, start, size):
            """Swaps a portion of data to the front of the array"""
            data[:size], data[start:start+size] = data[start:start+size], data[:size]
        log.info("%d pos and %d neg articles", len(self.pos), len(self.neg))
        pos = self.pos.copy()
        neg = self.neg.copy()
        if self.randomise:
            random.shuffle(pos)
            random.shuffle(neg)
        pstarts, psizes = self.partitionSizes(len(pos), self.nfold)
        nstarts, nsizes = self.partitionSizes(len(neg), self.nfold)
        pscores = numpy.zeros(len(pos), dtype=numpy.float32)
        nscores = numpy.zeros(len(neg), dtype=numpy.float32)
        pcounts = article.countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = article.countFeatures(self.numfeats, self.featdb, self.neg)
        for fold, (pstart,psize,nstart,nsize) in enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            log.debug("Carrying out fold number %d", fold)
            if statfile is not None:
                statfile.update(fold+1, self.nfold)
            moveToFront(pos, pstart, psize)
            moveToFront(neg, nstart, nsize)
            #print "POS : ", self.pos[:psize], self.pos[psize:]
            #print "NEG : ", self.neg[:nsize], self.neg[nsize:]
            temp_pcounts = pcounts - article.countFeatures(self.numfeats, self.featdb, self.pos[:psize])
            temp_ncounts = ncounts - article.countFeatures(self.numfeats, self.featdb, self.neg[:nsize])
            #print "PCNT: ", temp_pcounts
            #print "NCNT: ", temp_ncounts
            termscores, pfreqs, nfreqs = scoring.calculateFeatureScores(
                temp_pcounts, temp_ncounts, len(pos)-psize, len(neg)-nsize,
                self.pseudocount, self.daniel, self.mask)
            #print "TER: ", termscores
            #print
            pscores[pstart:pstart+psize] = [numpy.sum(termscores[self.featdb[d]]) for d in pos[:psize]]
            nscores[nstart:nstart+nsize] = [numpy.sum(termscores[self.featdb[d]]) for d in neg[:nsize]]
        return pscores, nscores
    
    def cross_validate_old(self, statfile=None):
        """Perform n-fold validation and return the raw performance measures"""
        def moveToFront(data, start, size):
            """Swaps a portion of data to the front of the array"""
            data[:size], data[start:start+size] = data[start:start+size], data[:size]
        log.info("%d pos and %d neg articles", len(self.pos), len(self.neg))
        pos = self.pos.copy()
        neg = self.neg.copy()
        if self.randomise:
            random.shuffle(pos)
            random.shuffle(neg)
        pstarts, psizes = self.partitionSizes(len(pos), self.nfold)
        nstarts, nsizes = self.partitionSizes(len(neg), self.nfold)
        pscores = numpy.zeros(len(pos), dtype=numpy.float32)
        nscores = numpy.zeros(len(neg), dtype=numpy.float32)
        for fold, (pstart,psize,nstart,nsize) in enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            log.debug("Carrying out fold number %d", fold)
            if statfile is not None:
                statfile.update(fold+1, self.nfold)
            moveToFront(pos, pstart, psize)
            moveToFront(neg, nstart, nsize)
            #print "POS: ", self.pos[:psize], self.pos[psize:]
            #print "NEG: ", self.neg[:nsize], self.neg[nsize:]
            pcounts = article.countFeatures(self.numfeats, self.featdb, self.pos[psize:])
            ncounts = article.countFeatures(self.numfeats, self.featdb, self.neg[nsize:])
            #print "PCNT: ", pcounts
            #print "NCNT: ", ncounts
            termscores, pfreqs, nfreqs = scoring.calculateFeatureScores(
                pcounts, ncounts, len(pos)-psize, len(neg)-nsize,
                self.pseudocount, self.daniel, self.mask)
            #print "TER: ", termscores
            #print
            pscores[pstart:pstart+psize] = [numpy.sum(termscores[self.featdb[d]]) for d in pos[:psize]]
            nscores[nstart:nstart+nsize] = [numpy.sum(termscores[self.featdb[d]]) for d in neg[:nsize]]
        return pscores, nscores

    def leave_out_one_validate(self, statfile=None):
        """Carries out leave-out-one validation, returning the resulting scores.

        @note: Does not support 'daniel' scoring method
        """
        pcounts = article.countFeatures(self.numfeats, self.featdb, self.pos)
        ncounts = article.countFeatures(self.numfeats, self.featdb, self.neg)
        #print "PCNT: ", pcounts
        #print "NCNT: ", ncounts
        pscores = numpy.zeros(len(self.pos), dtype=numpy.float32)
        nscores = numpy.zeros(len(self.neg), dtype=numpy.float32)
        pdocs = len(self.pos)
        ndocs = len(self.neg)
        ps = self.pseudocount
        for idx, pdoc in enumerate(self.pos):
            if statfile and idx % ((pdocs+ndocs)/20) == 0:
                statfile.update(idx, pdocs+ndocs)
            for f in self.featdb[pdoc]:
                if self.mask is None or self.mask[f] == False:
                    #print "TERM ", f, math.log(((pcounts[f]-1+ps)/(pdocs-1+2*ps))/((ncounts[f]+ps)/(ndocs+2*ps)))
                    pscores[idx] += math.log(((pcounts[f]-1+ps)/(pdocs-1+2*ps))/((ncounts[f]+ps)/(ndocs+2*ps)))
                    #pscores[idx] += math.log(((pcounts[f]+ps)/(pdocs+2*ps))/((ncounts[f]+ps)/(ndocs+2*ps)))
            #print pscores[idx]
        #print pscores
        for idx, ndoc in enumerate(self.neg):
            if statfile and (pdocs+idx) % ((pdocs+ndocs)/20) == 0:
                statfile.update(pdocs+idx, pdocs+ndocs)
            for f in self.featdb[ndoc]:
                if self.mask is None or self.mask[f] == False:
                    #print "TERM ", f, math.log(((pcounts[f]+ps)/(pdocs+2*ps))/((ncounts[f]-1+ps)/(ndocs-1+2*ps)))
                    nscores[idx] += math.log(((pcounts[f]+ps)/(pdocs+2*ps))/((ncounts[f]-1+ps)/(ndocs-1+2*ps)))
                    #nscores[idx] += math.log(((pcounts[f]+ps)/(pdocs+2*ps))/((ncounts[f]+ps)/(ndocs+2*ps)))
            #print nscores[idx]
        #print nscores
        return pscores, nscores

class PerformanceStats:
    """Calculates and stores performance statistics. There are arrays
    of TPR, FPR, PPV, FM, FMa for graphing, and counts P, N, A.

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

        @return: TPR, FPR, PPV, FM, FMa, ROC_area, PR_area, ThresholdIndex, Threshold, TP, FN, TN, FP, 
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
        idx, threshold =_.tuneThreshold()
        _.tunedStatistics(idx)

    def makeCountVectors(self):
        """Calculates arrays of TP, TN, FP, FN by iterating over pscores"""
        _ = self
        z = numpy.zeros(_.P, dtype=numpy.float32)
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
        z = numpy.zeros(_.P, dtype=numpy.float32)
        # TPR is recall
        _.TPR = _.TP / _.P
        # FPR is 1-specificity
        _.FPR = _.FP / _.N
        # PPV is precision
        _.PPV = _.TP / (_.TP + _.FP) 
        _.PPV[_.TP+_.FP == 0] = 1.0
        # FM is F-Measure
        prec = _.PPV
        rec = _.TPR
        _.FM = 2 * rec * prec / (rec + prec) 
        # FMa is the alpha-weighted F-Measures
        _.FMa = 1 / (_.alpha / prec + (1 - _.alpha) / rec)
        return _.TPR, _.FPR, _.PPV, _.FM, _.FMa

    def makeCurveAreas(self):
        """Calculate areas under ROC and precision-recall curves"""
        _ = self
        # Vectorised calculations
        import numpy as n
        _.ROC_area = n.sum(0.5 * (_.TPR[1:]+_.TPR[:-1]) * n.abs(n.diff(_.FPR))) + (1.0-max(_.FPR))
        _.PR_area = n.sum(_.PPV) / _.P
        # Non-vector calculation of ROC area
        ROC_area = 0
        for i in xrange(1,_.P):
            ROC_area += 0.5 * (_.TPR[i]+_.TPR[i-1]) * abs(_.FPR[i]-_.FPR[i-1])
        ROC_area += (1.0-max(_.FPR))
        # Non-vector calculation of PR area
        PR_area = 0
        for xi in xrange(0,_.P):
            PR_area += _.PPV[xi] / _.P
        print _.ROC_area, ROC_area
        print _.PR_area, PR_area
        return _.ROC_area, _.PR_area

    def tuneThreshold(self):
        """Return index into arrays which results in the best
        alpha-weighted F-Measure, and the score threshold which that
performance statistics        represents.
        """
        best_idx = 0
        for idx in xrange(self.P):
            if self.FMa[idx] > self.FMa[best_idx]:
                best_idx = idx
        return best_idx, self.pscores[best_idx]

    def tunedStatistics(self, index):
        """Object contains performance statistics for a particular tuned threshold.
    
        Instance variables include:
        index, threshold,
        TP, FP, TN, FN,
        P, N, A, T, F,
        TPR, FNR, TNR, FPR, PPV, NPV,
        accuracy, prevalence,
        recall, precision,
        fmeasure, fmeasure_alpha, fmeasure_max,
        enrichment
        """
        threshold = self.pscores[index]
        TP = self.TP[index]
        TN = self.TN[index]
        FP = self.FP[index]
        FN = self.FN[index]
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
        class Empty: pass
        self.tuned = Empty()
        self.tuned.__dict__.update(locals())
        del self.tuned.self
        del self.tuned.Empty
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
        daniel = c.dodaniel,
        exclude_types = c.exclude_types
        )
    feature_info.writeFeatureScoresCSV(codecs.open(rd/c.term_scores_name,"wb","utf-8"))
    #self.plotHistograms(gp, rd/c.hist_img, pscores, nscores, p.tuned.threshold)
    plotting.plotPDF(
        gp, rd/c.hist_img, p.pscores, p.nscores, p.tuned.threshold)
    plotting.plotROC(
        gp, rd/c.roc_img, p.FPR, p.TPR, p.tuned.FPR)
    plotting.plotPrecisionRecall(
        gp, rd/c.p_vs_r_img, p.TPR, p.PPV, p.tuned.TPR)
    plotting.plotPrecisionRecallFmeasure(
        gp, rd/c.pr_vs_score_img, p.pscores, p.TPR, p.PPV, p.FM, p.FMa, p.tuned.threshold)
    # Write main index file for validation output
    templates.validation.run(dict(
        p = p,
        t = p.tuned,
        f = feature_info,
        c = configuration,
        ), outputfile=file(rd/c.index_file, "w"))
    c.stylesheet.copy(rd/"style.css")
