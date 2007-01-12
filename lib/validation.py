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
import random
import time
import sys
#Local
from Gnuplot import Gnuplot, Data
import numpy
from path import path
#MScanner
import article
import scoring
import templates

class Validator:
    """Cross-validated performance statistics
    """

    def __init__(
        self,
        featmap,
        featdb,
        posids,
        negids,
        nfold,
        pseudocount=0.1,
        fm_tradeoff=0.0,
        dataset="",
        daniel=False,
        genedrug_articles=None,
        exclude_feats=None
        ):
        """Initialise validator
        @param featmap: Maping of feature id to (feature string, feature type)
        @param featdb: Mapping of doc id to list of feature ids
        @param posids: Array of positive doc ids
        @param negids: Array of negative doc ids
        @param nfold: Number of folds in cross-validation
        @param pseudocount, daniel: Parameters for scoring.getTermScores
        @param dataset: Title of the dataset being processed
        @param fm_tradeoff: Proportion of maximum F-Measure to trade in favour of precision
        @param genedrug_articles: Set of doc ids with gene-drug co-occurrences
        @param exclude_feats: Feature types to exclude for scoring purposes
        """
        self.featmap = featmap
        self.featdb = featdb
        self.pos = posids
        self.neg = negids
        self.nfold = nfold
        self.pseudocount = pseudocount
        self.fm_tradeoff = fm_tradeoff
        self.dataset = dataset
        self.daniel = daniel
        self.genedrug_articles = genedrug_articles
        self.exclude_feats = exclude_feats
        self.gnuplot = Gnuplot(debug=1)

    def articleIsPositive(self, docid, score, threshold):
        """Classifies an article as positive or negative based on score threshold"""
        if self.genedrug_articles is None or docid in self.genedrug_articles:
            return score >= threshold
        else:
            return False

    def validate(self, statfile=None):
        """Perform n-fold validation and return the raw performance measures"""
        def partitionSizes(nitems, nparts):
            """Returns start indeces and lengths of partitions"""
            base, rem = divmod(nitems, nparts)
            sizes = base * numpy.ones(nparts, dtype=numpy.int32)
            sizes[:rem] += 1
            starts = numpy.zeros(nparts, dtype=numpy.int32)
            starts[1:] = numpy.cumsum(sizes[:-1])
            return starts, sizes
        def moveToFront(data, start, size):
            """Swaps a portion of data to the front of the array"""
            data[:size], data[start:start+size] = data[start:start+size], data[:size]
        log.info("%d pos and %d neg articles", len(self.pos), len(self.neg))
        random.shuffle(self.pos)
        random.shuffle(self.neg)
        pstarts, psizes = partitionSizes(len(self.pos), self.nfold)
        nstarts, nsizes = partitionSizes(len(self.neg), self.nfold)
        pscores = numpy.zeros(len(self.pos), dtype=numpy.float32)
        nscores = numpy.zeros(len(self.neg), dtype=numpy.float32)
        for fold, (pstart,psize,nstart,nsize) in enumerate(zip(pstarts,psizes,nstarts,nsizes)):
            log.debug("Carrying out fold number %d", fold)
            if statfile is not None:
                statfile.update(fold+1, self.nfold)
            moveToFront(self.pos, pstart, psize)
            moveToFront(self.neg, nstart, nsize)
            pfreqs = article.countFeatures(len(self.featmap), self.featdb, self.pos[psize:])
            nfreqs = article.countFeatures(len(self.featmap), self.featdb, self.neg[nsize:])
            termscores = scoring.getTermScores(
                pfreqs, nfreqs, len(self.pos)-psize, len(self.neg)-nsize,
                self.pseudocount, self.daniel, self.featmap, self.exclude_feats)
            pscores[pstart:pstart+psize] = [scoring.scoreDocument(self.featdb[d], termscores) for d in self.pos[:psize]]
            nscores[nstart:nstart+nsize] = [scoring.scoreDocument(self.featdb[d], termscores) for d in self.neg[:nsize]]
        return pscores, nscores

    def plotPDF(self, fname, pdata, ndata, threshold):
        """Plot PDFs for pos/neg scores, with line to mark threshold""" 
        def kernelPDF(values, npoints=512):
            """Given 1D values, return an approximate probability density function
            @param values: Sorted list of floats representing the sample
            @param npoints: Number of equal-spaced points at which to estimate the PDF
            @return: (xvalues,yvalues) for y=f(x) of the pdf.
            """
            from scipy import stats
            points = numpy.linspace(values[0], values[-1], npoints)
            density = stats.kde.gaussian_kde(numpy.array(values)).evaluate(points)
            return points, density
        px, py = kernelPDF(pdata)
        nx, ny = kernelPDF(ndata)
        g = self.gnuplot
        g.reset()
        g.ylabel("Probability Density")
        g.xlabel("Article score")
        g.title("Score Densities")
        g("set terminal png")
        g("set output '%s'" % fname)
        threshold_height = max(chain(py, ny))
        g.plot(Data([threshold, threshold], [0, threshold_height], title="threshold", with="lines"),
               Data(px, py, title="Positives", with="lines"),
               Data(nx, ny, title="Negatives", with="lines"))
        
    def plotHistograms(self, fname, pdata, ndata, threshold):
        """Plot histograms for pos/neg scores, with line to mark threshold""" 
        py, px = numpy.histogram(pdata)
        ny, nx = numpy.histogram(ndata)
        g = self.gnuplot
        g.ylabel("Number of articles")
        g.xlabel("Score of article")
        g.title("Score Histograms")
        g("set terminal png")
        g("set output '%s'" % fname)
        threshold_height = max(chain(py, ny))
        g.plot(Data([threshold, threshold], [0, threshold_height], title="threshold", with="lines"),
               Data(px, py, title="Positives", with="histeps"),
               Data(nx, ny, title="Negatives", with="histeps"))

    @staticmethod
    def performanceStats(pscores, nscores, fm_tradeoff):
        """Derive arrays for TPR, FPR, PPV and FM over positive score
        points, and find peak performance

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @param fm_tradeoff: Proportion of maximum F-Measure to trade for increase precision.

        @note: Threshold is tuned to maximise precision subject to the
        F-Measure being at least 80% of the best F-Measure attained.

        @return: TPR, FPR, PPV, FM, ROC_area, PR_area, ThresholdIndex, Threshold, TP, FN, TN, FP, 
        """
        # Initialisation
        pscores.sort()
        nscores.sort()
        P, N = len(pscores), len(nscores)
        TPR = [ 0 for xi in xrange(P) ] # recall
        FPR = [ 0 for xi in xrange(P) ] # 1-specificity
        PPV = [ 0 for xi in xrange(P) ] # precision
        FM  = [ 0 for xi in xrange(P) ] # F-measure
        TN  = 1                         # True negatives
        maxFM_xi = 0                    # Maximum of F-measure
        best_xi, best_TP, best_FN, best_TN, best_FP = 0, P, 0, N, 0
        ROC_area, PR_area = 0, 0
        # Calculate stats for each choice of threshold
        for xi in xrange(P):
            threshold = pscores[xi]
            while (nscores[TN-1] < threshold) and (TN < N):
                TN += 1
            # xi+1 is how many positives we called negative
            FN = xi+1
            # TP+FN = P
            TP = P - FN
            # TN+FP = N
            FP = N - TN
            # TPR = TP/P
            TPR[xi] = TP/P
            # FPR = FP/N = 1 - TN/N = 1 - specificity
            FPR[xi] = FP/N 
            # PPV = TP/(TP+FP) = precision
            PPV[xi] = 0
            if TP+FP > 0:
                PPV[xi] = TP/(TP+FP) 
            # F-Measure = 2*recall*precision/(recall+precision)
            FM[xi] = 0
            if PPV[xi] > 0 and TPR[xi] > 0:
                FM[xi] = 2*TPR[xi]*PPV[xi]/(PPV[xi]+TPR[xi])
            if xi > 0:
                # Use trapezoidal rule to integrate ROC and PR curves
                ROC_area += 0.5*abs(TPR[xi]+TPR[xi-1])*abs(FPR[xi]-FPR[xi-1])
                PR_area += 0.5*abs(PPV[xi]+PPV[xi-1])*abs(TPR[xi]-TPR[xi-1])
                #print "PPV=%g, TPR=%g, FPR=%g, ROC=%g, PR=%g" % (PPV[xi], TPR[xi], FPR[xi], ROC_area, PR_area)
            # Track maximum F-Measure
            if FM[xi] > FM[maxFM_xi]:
                maxFM_xi = xi
            # Tune to maximise PPV subject to at least 90% of maximum F-Measure
            if PPV[xi] > PPV[best_xi] and FM[xi] >= FM[maxFM_xi]*(1-fm_tradeoff):
                best_xi, best_TP, best_FN, best_TN, best_FP = xi, TP, FN, TN, FP
            #print "thresh = %g, TPR = %d/%d = %.1e, FPR = %d/%d = %.1e" % (threshold, TP, P, TP/P, FP, N, FP/N)
        ROC_area += (1-max(FPR))
        # Return tuned results
        return TPR, FPR, PPV, FM, ROC_area, PR_area, best_xi, pscores[best_xi], best_TP, best_FN, best_TN, best_FP

    def plotROC(self, roc, TPR, FPR, thresh_idx):
        """ROC curve (TPR vs FPR)

        @param roc: Path to output file for ROC curve
        """
        g = self.gnuplot
        g.reset()
        g.ylabel("True Positive Rate (TPR)")
        g.xlabel("False Positive Rate (FPR)")
        g.title("ROC curve (TPR vs FPR)")
        g("set terminal png")
        g("set output '%s'" % roc)
        g.plot(Data(FPR, TPR, title="TPR", with="lines"),
               Data([FPR[thresh_idx], FPR[thresh_idx]], [0,1.0], title="threshold", with="lines"))

    def plotPR(self, p_vs_r, TPR, PPV, thresh_idx):
        """Precision vs recall graph

        @param p_vs_r: Path to output file for precision-recall curve
        """
        g = self.gnuplot
        g.reset()
        g.ylabel("Precision")
        g.xlabel("Recall")
        g.title("Precision vs Recall")
        g("set terminal png")
        g("set output '%s'" % p_vs_r)
        g.plot(Data(TPR, PPV, title="Precision", with="lines", smooth="csplines"),
               Data([TPR[thresh_idx], TPR[thresh_idx]], [0,1.0], title="threshold", with="lines"))

    def plotPRF(self, pr_vs_score, pscores, TPR, PPV, FM, threshold):
        """Precision, Recall, F-Measure vs threshold graph

        @parav pr_vs_score: Path to output file for precision,recall vs threshold
        """
        g = self.gnuplot
        g.reset()
        g.ylabel("Precision, Recall, F-Measure")
        g.xlabel("Threshold Score")
        g.title("Precision and Recall vs Threshold")
        g("set terminal png")
        g("set output '%s'" % pr_vs_score)
        g.plot(Data(pscores, TPR, title="Recall", with="lines"),
               Data(pscores, PPV, title="Precision", with="lines"),
               Data(pscores, FM, title="F-Measure", with="lines"),
               Data([threshold, threshold], [0,1], title="threshold", with="lines"))

    def report(self, pscores, nscores, prefix, stylesheet):
        """Write a full validation report

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @param prefix: Directory for validation report validation report

        @param stylesheet: Path to CSS stylesheet for the report templates
        """
        # Output files
        terms_csv = "termscores.csv"
        hist_img = "histogram.png"
        roc_img = "roc.png"
        p_vs_r_img = "prcurve.png"
        pr_vs_score_img = "prscore.png"
        mainfile = "index.html"
        graph_pickle = "graphs.pickle"
        # Performance tuning, save arrays for graphing
        aTPR, aFPR, aPPV, aFM, ROC_area, PR_area, thresh_idx, threshold, TP, FN, TN, FP  = self.performanceStats(pscores, nscores, self.fm_tradeoff)
        cPickle.dump(dict(TPR=aTPR,FPR=aFPR,PPV=aPPV,FM=aFM), file(graph_pickle, "wb"))
        #self.plotHistograms(hist_img, pscores, nscores, threshold)
        self.plotPDF(prefix/hist_img, pscores, nscores, threshold)
        self.plotROC(prefix/roc_img, aTPR, aFPR, thresh_idx)
        self.plotPR(prefix/p_vs_r_img, aTPR, aPPV, thresh_idx)
        self.plotPRF(prefix/pr_vs_score_img, pscores, aTPR, aPPV, aFM, threshold)
        # Calculate and write feature scores to CSV
        pfreqs = article.countFeatures(len(self.featmap), self.featdb, self.pos)
        nfreqs = article.countFeatures(len(self.featmap), self.featdb, self.neg)
        pdocs, ndocs = len(self.pos), len(self.neg)
        termscores = scoring.getTermScores(
            pfreqs, nfreqs, pdocs, ndocs, self.pseudocount, 
            self.daniel, self.featmap, self.exclude_feats)
        scoring.writeFeatureScoresCSV(codecs.open(prefix/terms_csv,"wb","utf-8"), self.featmap, termscores)
        feature_stats = scoring.featureStatistics(self.featmap, termscores, pdocs, ndocs)
        # Calculate performance measures
        P = TP+FN
        N = TN+FP
        A = TP+FN+TN+FP # A = P + N = T + F
        T = TP+TN
        F = FP+FN
        TPR, FNR, TNR, FPR, PPV, NPV = 0, 0, 0, 0, 0, 0
        accuracy, prevalence, enrichment, fmeasure = 0, 0, 0, 0
        if TP+FN != 0:
            TPR = TP/(TP+FN) # TPR = TP/P = sensitivity = recall
            FNR = FN/(TP+FN) # FNR = FN/P = 1 - TP/P = 1-sensitivity = 1-recall
        if TN+FP != 0:
            TNR = TN/(TN+FP) # TNR = TN/N = specificity
            FPR = FP/(TN+FP) # FPR = FP/N = 1 - TN/N = 1-specificity
        if TP+FP != 0:
            PPV = TP/(TP+FP) # PPV = precision
        if TN+FN != 0:
            NPV = TN/(TN+FN) # NPV
        if A != 0:
            accuracy = (TP+TN)/A   # acc  = T/A
            prevalence = (TP+FN)/A # prev = P/A
        recall = TPR
        precision = PPV
        if recall > 0 and precision > 0:
            fmeasure = 2*recall*precision/(recall+precision)
        if prevalence > 0:
            enrichment = precision / prevalence
        # Write main index file for validation output
        dataset = self.dataset
        nfold = self.nfold
        templates.validation.run(locals(), outputfile=file(prefix/mainfile, "w"))
        stylesheet.copy(prefix / "style.css")        

