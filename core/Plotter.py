"""Plotting functions for all graphs produced in cross validation."""

from __future__ import division
from Gnuplot import Data, Gnuplot
import logging as log
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


class Plotter(Gnuplot):
    """Implements the plots used in MScanner
    
    When adding a new analysis, the plotting function for its output
    graphs should be added here.
    
    All functions take the fname parameter, which is the path to 
    the PNG file to which the graph will be written.
    """


    def plot_predictions(g, fname, predicted_low, predicted_high):
        """Given L{PredictedMetrics} instance, plot the predicted query
        performance (TPR and PPV vs number of results).
        
        @note: Two predictions are given, corresponding to upper and lower
        bound guesses at the number of 
        """
        log.debug("Plotting prediction curve to %s", fname)
        g.reset()
        g.title("Prediction curve (TPR, PPV vs # results)")
        g.ylabel("TPR, PPV")
        g.xlabel("Number of results")
        g("set terminal png")
        g("set output '%s'" % fname)
        # Only plot until 95% recall in both cases
        pL = predicted_low
        pH = predicted_high
        seg = (pL.TPR < 0.95) | (pH.TPR < 0.95) 
        g.plot(
            Data(pL.results[seg], pL.TPR[seg], title="TPR low", with="lines"),
            Data(pL.results[seg], pL.PPV[seg], title="PPV low", with="lines"),
            Data(pH.results[seg], pH.TPR[seg], title="TPR high", with="lines"),
            Data(pH.results[seg], pH.PPV[seg], title="PPV high", with="lines"))

    
    def plot_roc(g, fname, FPR, TPR, marker_FPR):
        """ROC curve (TPR vs FPR)"""
        log.debug("Plotting ROC curve to %s", fname)
        g.reset()
        g.title("ROC curve (TPR vs FPR)")
        g.ylabel("True Positive Rate (TPR)")
        g.xlabel("False Positive Rate (FPR)")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(FPR, TPR, title="TPR", with="lines"),
               Data([marker_FPR, marker_FPR], [0,0.99], title="threshold", with="lines"))


    def plot_precision(g, fname, TPR, PPV, marker_TPR):
        """Precision vs recall"""
        log.debug("Plotting Precision-Recall curve to %s", fname)
        g.reset()
        g.title("Precision vs Recall")
        g.ylabel("Precision")
        g.xlabel("Recall")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(TPR, PPV, title="Precision", with="lines", smooth="csplines"),
               Data([marker_TPR, marker_TPR], [0,0.99], title="threshold", with="lines"))


    def plot_fmeasure(g, fname, pscores, TPR, PPV, FM, FMa, threshold):
        """Precision, Recall, F-Measure vs threshold"""
        log.debug("Plotting F-Measure curve to %s", fname)
        g.reset()
        g.title("Precision and Recall vs Threshold")
        g.ylabel("Precision, Recall, F-Measure, F-Measure Alpha")
        g.xlabel("Threshold Score")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(pscores, TPR, title="Recall", with="lines"),
               Data(pscores, PPV, title="Precision", with="lines"),
               Data(pscores, FM, title="F1 Measure", with="lines"),
               Data(pscores, FMa, title="F Measure", with="lines"),
               Data([threshold, threshold], [0,0.99], title="threshold", with="lines"))


    @staticmethod
    def bincount(data):
        """Calculate the best number of histogram bins for the data
        
        Uses the formula M{K = R/(2*IQR*N^(-1/3))}
        
        @param data: Array of numbers, sorted in increasing order.
        """
        N = len(data) # Number of data points
        IQR = data[3*N//4] - data[N//4] # Inter-Quartile Range
        R = data[-1] - data[0] # Range
        bins = R//(2*IQR*N**(-1/3)) # Number of bins
        #print N, IQR, R, bins
        return min(150, max(10, bins))


    def plot_score_histogram(g, fname, pdata, ndata, threshold):
        """Histograms for pos and neg scores, with line to mark threshold""" 
        log.debug("Plotting article score histogram to %s", fname)
        from itertools import chain
        py, px = nx.histogram(pdata, bins=g.bincount(pdata), normed=True)
        zy, zx = nx.histogram(ndata, bins=g.bincount(ndata), normed=True)
        g.reset()
        g("set terminal png")
        g("set output '%s'" % fname)
        g.title("Score Histograms")
        g.xlabel("Article Score")
        g.ylabel("Histogram Mass")
        ## Commented out arrow - rather adding a real line
        #g("set arrow from %f,0 to %f,%f nohead lw 4 " % (
        #    threshold, threshold, max(chain(py,ny))))
        g("set style fill solid 1.0")
        threshold_height = max(chain(py, zy))
        g.plot(Data(px, py, title="Positives", with="boxes"),
               Data(zx, zy, title="Negatives", with="boxes"),
               Data([threshold, threshold], [0, threshold_height], 
                    title="threshold", with="lines lw 3"))


    def plot_feature_histogram(g, fname, scores):
        """Histogram for feature scores
        
        @param scores: List with scores of each feature"""
        log.debug("Plotting feature score histogram to %s", fname)
        sscores = scores.copy()
        sscores.sort()
        y, x = nx.histogram(scores, bins=g.bincount(sscores))
        g.reset()
        g.title("Feature Score Histogram")
        g.xlabel("Feature Score")
        g.ylabel("Histogram Mass")
        g("set logscale y")
        g("set terminal png")
        g("set output '%s'" % fname)
        g("set style fill solid 1.0")
        g.plot(Data(x, y, with="boxes"))



class RetrievalPlot(Gnuplot):
    """Plotting function for comparative retrieval test"""
    
    def plot_retrieved_positives(g, fname, nretrieved, total):
        """Proportion of testing set retrieved, versus result rank"""
        log.debug("Plotting Retrieval curve to %s", fname)
        g.reset()
        g.title("Retrieval Test")
        g.ylabel("Proportion Retrieved (recall)")
        g.xlabel("Rank")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(range(0,len(nretrieved)), nretrieved / total, with="lines"))



def DensityPlotter(Plotter):
    """Adds plotting of estimated Probability Density Functions for
    article and feature scores.
    
    @deprecated: These methods are too computationally
    expensive to use interactively (about 40 seconds per graph).
    """


    @staticmethod
    def gaussian_kernel_pdf(values, npoints=512):
        """Given 1D values, return the probability density function
        
        @param values: Sorted list of floats representing the sample
        
        @param npoints: Number of equal-spaced points at which to estimate the PDF
        
        @return: (xvalues, yvalues) for y=f(x) of the pdf.
        """
        from scipy import stats
        points = nx.linspace(values[0], values[-1], npoints)
        density = stats.kde.gaussian_kde(nx.array(values)).evaluate(points)
        return points, density


    def plot_score_density(g, fname, pdata, ndata, threshold):
        """Probability density of pos and neg scores, with line to mark threshold
    
        @param pdata: Scores of positive documents
        @param ndata: Scores of negative documents
        @param threshold: Threshold score for counting a document positive
        """ 
        from itertools import chain
        log.debug("Plotting article score density to %s", fname)
        px, py = g.gaussian_kernel_pdf(pdata)
        zx, zy = g.gaussian_kernel_pdf(ndata)
        overlap = calculateOverlap(px, py, zx, zy)
        g.reset()
        g.title("Article Score Densities")
        g.ylabel("Probability Density")
        g.xlabel("Article score")
        g("set terminal png")
        g("set output '%s'" % fname)
        threshold_height = max(chain(py, zy))
        g.plot(Data([threshold, threshold], [0, threshold_height], 
                    title="threshold", with="lines"),
               Data(px, py, title="Positives", with="lines"),
               Data(zx, zy, title="Negatives", with="lines"))
        return overlap


    def plot_feature_density(g, fname, scores):
        """Probability density function for feature scores"""
        log.debug("Plotting feature score density to %s", fname)
        x, y = g.gaussian_kernel_pdf(scores, npoints=1024)
        g.reset()
        g.title("Feature Score Density")
        g.xlabel("Feature Score")
        g.ylabel("Probability Density")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(x, y, with="lines"))


