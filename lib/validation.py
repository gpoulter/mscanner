"""Cross-validation and performance statistics module

@author: Graham Poulter
                                   

Validator -- Perform cross-validation and output performance statistics

"""

from __future__ import division
from article import TermCounts
from Gnuplot import Gnuplot, Data
from heapq import heapreplace
from itertools import chain
import logging as log
from path import path
from random import randint, seed
import templates
import time
from scoring import getTermScores, scoreDocument, writeTermScoresCSV, writeTermScoresHTML
import sys
import warnings
warnings.filterwarnings(action='ignore', category=FutureWarning)

class Validator:
    """Cross-validated performance statistics

    Trains on 90% of positives and negatives, and use the other 10%
    for testing.
    """

    def __init__(
        self,
        meshdb,
        featdb,
        posids,
        negids,
        nfold,
        pseudocount=0.1,
        daniel=False,
        genedrug_articles=None ):
        """Initialise validator
        @param meshdb: Maping of term id to term name
        @param featdb: Maping of doc id to list of term ids
        @param posids: Set of positive doc ids
        @param negids: Set of negative doc ids
        @param nfold: Number of folds in cross-validation.
        @param pseudocount, daniel: Passed to Scoring.getTermScores
        @param genedrug_articles: Set of doc ids with gene-drug co-occurrences
        """
        self.meshdb = meshdb
        self.featdb = featdb
        self.pos = posids
        self.neg = negids
        self.nfold = nfold
        self.pseudocount = pseudocount
        self.daniel = daniel
        self.genedrug_articles = genedrug_articles
        self.gnuplot = Gnuplot(debug=1)

    def articleIsPositive( self, docid, score, threshold ):
        """Classifies an article as positive or negative based on score threshold"""
        if self.genedrug_articles is None or docid in self.genedrug_articles:
            return score >= threshold
        else:
            return False

    @staticmethod
    def partition(items, nfold):
        """Partitions indeces into data items for cross-validation

        @param items: Set of objects to be partitioned
        @param nfold: Number of partitions to use
        @return: List of sets, each set being length len(items)/nfold
        """
        tests = [ set() for n in range(nfold) ]
        ct = 0                  # current test
        available = list(items) # items to assign
        N = len(available)      # number of items
        i = 0                   # current iteration
        size = N                # number of items left
        while i < N:
            r = randint(0,size-1)            # choose a random item
            tests[ct].add(available[r])      # add item to current set
            available[r] = available[size-1] # over write with last item
            i += 1                           # next iteration
            ct = (ct+1) % len(tests)         # next test
            size = size - 1                  # virtual deletion
        return tests

    def validate(self):
        """Perform n-fold validation and return the raw performance measures"""
        positives = self.pos
        negatives = self.neg
        log.info( "%d pos and %d neg articles", len(positives), len(negatives) )
        ptests = self.partition(positives, self.nfold)
        ntests = self.partition(negatives, self.nfold)
        pscores, nscores, threshold = [],[],0
        for fold, (ptest, ntest) in enumerate(zip(ptests, ntests)):
            log.debug("Carrying out fold number %d", fold)
            ptrain = positives - ptest
            ntrain = negatives - ntest
            pfreqs = TermCounts(self.featdb[d] for d in ptrain)
            nfreqs = TermCounts(self.featdb[d] for d in ntrain)
            termscores = getTermScores(pfreqs, nfreqs, self.pseudocount, self.daniel)
            pscores.extend(scoreDocument(self.featdb[d], termscores) for d in ptest)
            nscores.extend(scoreDocument(self.featdb[d], termscores) for d in ntest)
        return pscores, nscores

    @staticmethod
    def makeHistogram(data):
        """Returns a tuple (centers,freqs) of normalised frequencies
        at centers, given data points in data and a set bin width."""
        small = min(data)
        big = max(data)
        npoints = len(data)
        bin_width = min(npoints, (big-small)/20)
        centers = [ 0.0 for x in range(0,int((big+1e-5-small)/bin_width)) ]
        if len(centers) == 0:
            return [0], [0]
        for i in range(len(centers)):
            centers[i] = small+bin_width*(2*i+1)/2
        freqs = [ 0.0 for x in centers ]
        for x in data:
            pos = int((x-small)//bin_width)
            if pos == len(freqs): pos = len(freqs)-1
            freqs[pos] += 1/npoints
        return centers, freqs

    @staticmethod
    def kernelPDF(values, npoints=512):
        """Given 1D values, return an approximate probability density function

        @note: Uses SciPy, which uses Gaussian kernel

        @param values: Sorted list of floats representing the sample

        @param npoints: Number of equal-spaced points at which to estimate the PDF

        @return: (xvalues,yvalues) for y=f(x) of the pdf.
        """
        from scipy import stats
        import numpy
        points = numpy.linspace(values[0], values[-1], npoints)
        return points, stats.kde.gaussian_kde(numpy.array(values)).evaluate(points)

    def plotPDF(self, name, pdata, ndata, threshold):
        """Plot PDFs for pos/neg scores, with line to mark threshold""" 
        px, py = self.kernelPDF(pdata)
        nx, ny = self.kernelPDF(ndata)
        g = self.gnuplot
        g.reset()
        g.ylabel("Probability Density")
        g.xlabel("Article score")
        g.title("Score Densities")
        g("set terminal png")
        g("set output '%s'" % name)
        threshold_height = max(chain(py, ny))
        g.plot(Data([threshold, threshold], [0, threshold_height], title="threshold", with="lines"),
               Data(px, py, title="Positives", with="lines"),
               Data(nx, ny, title="Negatives", with="lines"))
        
    def plotHistograms(self, name, pdata, ndata, threshold):
        """Plot histograms for pos/neg scores, with line to mark threshold""" 
        px, py = self.makeHistogram(pdata)
        nx, ny = self.makeHistogram(ndata)
        g = self.gnuplot
        g.ylabel("Probability Histogram")
        g.xlabel("Article score")
        g.title("Score Histograms")
        g("set terminal png")
        g("set output '%s'" % name)
        threshold_height = max(chain(py, ny))
        g.plot(Data([threshold, threshold], [0, threshold_height], title="threshold", with="lines"),
               Data(px, py, title="Positives", with="histeps"),
               Data(nx, ny, title="Negatives", with="histeps"))

    def plotCurves(self, roc, p_vs_r, pr_vs_score, pscores, nscores):
        """Plot curves for ROC, precision-recall, and stats vs threshold

        @param roc: Path to output file for ROC curve

        @param p_vs_r: Path to output file for precision-recall curve

        @parav pr_vs_score: Path to output file for precision,recall vs threshold

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @return: Threshold, TP, FN, TN, FP, such that threshold is
        tuned to maximise precision subject to the F-Measure being at
        least 80% of the best F-Measure attained.
        """
        # Initialisation
        pscores.sort()
        nscores.sort()
        P = len(pscores)
        N = len(nscores)
        TPR = [ 0 for xi in xrange(P) ] # recall
        FPR = [ 0 for xi in xrange(P) ] # 1-specificity
        PPV = [ 0 for xi in xrange(P) ] # precision
        FM  = [ 0 for xi in xrange(P) ] # F-measure
        TN  = 1                         # True negatives
        maxFM_xi = 0                    # Maximum of F-measure
        best_xi = 0
        best_TP = 0
        best_FN = 0
        best_TN = 0
        best_FP = 0
        ROC_area = 0
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
            if xi > 0:
                # Use trapezoidal rule to integrate ROC curve
                ROC_area += 0.5*(TPR[xi]+TPR[xi-1])*(FPR[xi-1]-FPR[xi])
            # PPV = TP/(TP+FP) = precision
            PPV[xi] = 0
            if TP+FP > 0:
                PPV[xi] = TP/(TP+FP) 
            # F-Measure = 2*recall*precision/(recall+precision)
            FM[xi] = 0
            if PPV[xi] > 0 and TPR[xi] > 0:
                FM[xi] = 2*TPR[xi]*PPV[xi]/(PPV[xi]+TPR[xi])
            # Track maximum F-Measure
            if FM[xi] > FM[maxFM_xi]:
                maxFM_xi = xi
            # Tune to maximise PPV subject to at least 90% of maximum F-Measure
            if PPV[xi] > PPV[best_xi] and FM[xi] > 0.9*FM[maxFM_xi]:
                best_xi, best_TP, best_FN, best_TN, best_FP = xi, TP, FN, TN, FP
            #print "thresh = %g, TPR = %d/%d = %.1e, FPR = %d/%d = %.1e" % (threshold, TP, P, TP/P, FP, N, FP/N)
        g = self.gnuplot
        # Precision vs recall graph
        g.reset()
        g.ylabel("Precision")
        g.xlabel("Recall")
        g.title("Precision vs Recall")
        g("set terminal png")
        g("set output '%s'" % p_vs_r)
        g.plot(Data(TPR, PPV, title="Precision", with="lines"),
               Data([TPR[best_xi], TPR[best_xi]], [0,1.0], title="threshold", with="lines"))
        # ROC curve (TPR vs FPR)
        g.reset()
        g.ylabel("True Positive Rate (TPR)")
        g.xlabel("False Positive Rate (FPR)")
        g.title("ROC curve (TPR vs FPR)")
        g("set terminal png")
        g("set output '%s'" % roc)
        g.plot(Data(FPR, TPR, title="TPR", with="lines"),
               Data([FPR[best_xi], FPR[best_xi]], [0,1.0], title="threshold", with="lines"))
        # Precision, Recall, F-Measure vs threshold graph
        g.reset()
        g.ylabel("Precision, Recall, F-Measure")
        g.xlabel("Threshold Score")
        g.title("Precision and Recall vs Threshold")
        g("set terminal png")
        g("set output '%s'" % pr_vs_score)
        g.plot(Data(pscores, TPR, title="Recall", with="lines"),
               Data(pscores, PPV, title="Precision", with="lines"),
               Data(pscores, FM, title="F-Measure", with="lines"),
               Data([pscores[best_xi], pscores[best_xi]], [0,1], title="threshold", with="lines"))
        # Return tuned results
        return pscores[best_xi], best_TP, best_FN, best_TN, best_FP, ROC_area

    def publicationGraphs(self, pscores, nscores, prefix):
        """Draws graphs for publication

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @param prefix: Directory for validation report validation report
        """
        hist_img = prefix/"histogram_pub.png"
        roc_img = prefix/"roc_pub.png"
        p_vs_r_img = prefix/"prcurve_pub.png"
        threshold, TP, FN, TN, FP, ROC_area = self.plotCurves(roc_img, p_vs_r_img, pr_vs_score_img, pscores, nscores)
        self.plotPDF(hist_img, pscores, nscores, threshold)

    def report(self, pscores, nscores, prefix, stylesheet):
        """Write a full validation report

        @param pscores: Scores of positive articles

        @param nscores: Scores of negative articles

        @param prefix: Directory for validation report validation report

        @param stylesheet: Path to CSS stylesheet for the report templates
        """
        # Set up output files
        terms_csv = prefix/"termscores.csv"
        terms_html = prefix/"termscores.html"
        hist_img = prefix/"histogram.png"
        roc_img = prefix/"roc.png"
        p_vs_r_img = prefix/"prcurve.png"
        pr_vs_score_img = prefix/"prscore.png"
        mainfile = prefix/"index.html"

        # Graphs and performance tuning
        threshold, TP, FN, TN, FP, ROC_area = self.plotCurves(roc_img, p_vs_r_img, pr_vs_score_img, pscores, nscores)
        self.plotHistograms(hist_img, pscores, nscores, threshold)

        # Output term scores
        pfreqs = TermCounts(self.featdb[d] for d in self.pos)
        nfreqs = TermCounts(self.featdb[d] for d in self.neg)
        termscores = getTermScores( pfreqs, nfreqs, self.pseudocount, self.daniel)
        writeTermScoresCSV(file(terms_csv,"w"), self.meshdb, termscores, pfreqs, nfreqs)
        writeTermScoresHTML(file(terms_html,"w"), self.meshdb, termscores, pfreqs, nfreqs, self.pseudocount)

        # Calculate performance measures
        P = TP+FN
        N = TN+FP
        A = TP+FN+TN+FP # = P + N = T + F
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

        # Write main index file
        templates.validation.run(dict(
            TP=TP, TN=TN, FP=FP, FN=FN, P=P, N=N, A=A, T=T, F=F,
            TPR=TPR, FNR=FNR, TNR=TNR, FPR=FPR, PPV=PPV, NPV=NPV,
            ROC_area = ROC_area,
            accuracy = accuracy,
            prevalence = prevalence,
            enrichment = enrichment,
            recall = recall,
            threshold = threshold,
            precision = precision,
            fmeasure = fmeasure,
            terms_csv = terms_csv.basename(),
            terms_html = terms_html.basename(),
            hist_img = hist_img.basename(),
            roc_img = roc_img.basename(),
            p_vs_r_img = p_vs_r_img.basename(),
            pr_vs_score_img = pr_vs_score_img.basename(),
            ), outputfile=file(mainfile, "w"))
        stylesheet.copy(prefix / "style.css")
