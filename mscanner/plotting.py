"""Plotting functions for all graphs produced in cross validation."""

from __future__ import division

                                     
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

from Gnuplot import Data, Gnuplot
import logging as log
import numpy as n


def bincount(data):
    """Return best number of histogram bins for the data 
    
    @note: Uses the formula M{K = R/(2*IQR*N^(-1/3))}

    @note: Expects data to be sorted increasing order.
    """
    N = len(data) # Number of data points
    IQR = data[3*N//4] - data[N//4] # Inter-Quartile Range
    R = data[-1] - data[0] # Range
    bins = R//(2*IQR*N**(-1/3)) # Number of bins
    #print N, IQR, R, bins
    return min(150,max(10, bins))


def calculateOverlap(px, py, nx, ny):
    """Calculate overlap between two bell curves (curve p must be
    to the right of curve n)
    
    Procedure is first to interpolate both curves onto a common X-axis ranging
    from min(nx) to max(px), then find highest where pY and nY intersect,
    then create a curve using pY up to intersection, and nY thereafter. Then
    calculate the area underneath using trapezoidal rule."""
    from scipy.interpolate import interp1d
    from scipy.integrate import trapz
    X = n.linspace(n.min(nx), n.max(px), 1000)
    p_interp = interp1d(px, py, bounds_error=False, fill_value=0.0)
    n_interp = interp1d(nx, ny, bounds_error=False, fill_value=0.0)
    pY = p_interp(X)
    nY = n_interp(X)
    # Attempt to find point of intersection, 
    # but leave out interpolated sections with zero density
    diffs = n.absolute(pY-nY)
    diffs[pY==0] = 1
    diffs[nY==0] = 1
    interidx = n.nonzero(diffs == n.min(diffs))[0][0]
    iY = n.concatenate((pY[:interidx],nY[interidx:]))
    area = trapz(iY,X)
    include = iY != 0
    ingraph = False
    for idx in xrange(len(iY)):
        if not ingraph and include[idx]:
            include[idx-1] = True
            ingraph = True
        if ingraph and not include[idx]:
            include[idx] = True
            break
    return area, X[include], iY[include]


def kernelPDF(values, npoints=512):
    """Given 1D values, return an approximate probability density function
    
    @param values: Sorted list of floats representing the sample
    
    @param npoints: Number of equal-spaced points at which to estimate the PDF
    
    @return: (xvalues,yvalues) for y=f(x) of the pdf.
    """
    from scipy import stats
    points = n.linspace(values[0], values[-1], npoints)
    density = stats.kde.gaussian_kde(n.array(values)).evaluate(points)
    return points, density



class Plotter(Gnuplot):
    
    def plotArticleScoreDensity(g, fname, pdata, ndata, threshold):
        """Probability density of pos and neg scores, with line to mark threshold
    
        @param fname: Filename to plot to
        @param pdata: Scores of positive documents
        @param ndata: Scores of negative documents
        @param threshold: Threshold score for counting a document positive
        """ 
        from itertools import chain
        log.debug("Plotting article score density to %s", fname)
        px, py = kernelPDF(pdata)
        nx, ny = kernelPDF(ndata)
        overlap = calculateOverlap(px, py, nx, ny)
        g.reset()
        g.title("Article Score Densities")
        g.ylabel("Probability Density")
        g.xlabel("Article score")
        g("set terminal png")
        g("set output '%s'" % fname)
        threshold_height = max(chain(py, ny))
        g.plot(Data([threshold, threshold], [0, threshold_height], title="threshold", with="lines"),
               Data(px, py, title="Positives", with="lines"),
               Data(nx, ny, title="Negatives", with="lines"))
        return overlap


    def plotFeatureScoreDensity(g, fname, scores):
        """Probability density function for feature scores"""
        log.debug("Plotting feature score density to %s", fname)
        x, y = kernelPDF(scores, npoints=1024)
        g.reset()
        g.title("Feature Score Density")
        g.xlabel("Feature Score")
        g.ylabel("Probability Density")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(x, y, with="lines"))


    def plotArticleScoreHistogram(g, fname, pdata, ndata, threshold):
        """Histograms for pos and neg scores, with line to mark threshold""" 
        log.debug("Plotting article score histogram to %s", fname)
        from itertools import chain
        py, px = n.histogram(pdata, bins=bincount(pdata), normed=True)
        ny, nx = n.histogram(ndata, bins=bincount(ndata), normed=True)
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
        threshold_height = max(chain(py, ny))
        g.plot(Data(px, py, title="Positives", with="boxes"),
               Data(nx, ny, title="Negatives", with="boxes"),
               Data([threshold, threshold], [0, threshold_height], 
                    title="threshold", with="lines lw 3"))


    def plotFeatureScoreHistogram(g, fname, scores):
        """Histogram for feature scores
        
        @param scores: List with scores of each feature"""
        log.debug("Plotting feature score histogram to %s", fname)
        sscores = scores.copy()
        sscores.sort()
        y, x = n.histogram(scores, bins=bincount(sscores))
        g.reset()
        g.title("Feature Score Histogram")
        g.xlabel("Feature Score")
        g.ylabel("Histogram Mass")
        g("set logscale y")
        g("set terminal png")
        g("set output '%s'" % fname)
        g("set style fill solid 1.0")
        g.plot(Data(x, y, with="boxes"))


    def plotROC(g, fname, FPR, TPR, marker_FPR):
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


    def plotPrecisionRecall(g, fname, TPR, PPV, marker_TPR):
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


    def plotPrecisionRecallFmeasure(g, fname, pscores, TPR, PPV, FM, FMa, threshold):
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


    def plotRetrievalGraph(g, fname, nretrieved, total):
        """Proportion of testing set retrieved, versus result rank"""
        log.debug("Plotting Retrieval curve to %s", fname)
        g.reset()
        g.title("Retrieval Test")
        g.ylabel("Proportion Retrieved (recall)")
        g.xlabel("Rank")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(range(0,len(nretrieved)), nretrieved / total, with="lines"))
