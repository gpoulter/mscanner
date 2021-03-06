"""Plotting functions for all graphs produced in cross validation."""

from __future__ import division
from Gnuplot import Data, Gnuplot
import logging
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


class Plotter:
    """Implements the plots used in MScanner
    
    When adding a new analysis, the plotting function for its output
    graphs should be added here.
    
    @note: All methods take an fname parameter, which is the path to the PNG
    file to which the graph will be written.
    
    @ivar overwrite: If False, we no-op rather than overwrite an already
    existing graph.
    
    @ivar gnuplot: The captive Gnuplot instance.
    """

    def __init__(self, overwrite=True):
        self.overwrite = overwrite
        self.gnuplot = Gnuplot()


    def interpolate(self, x, y):
        """Return 8000-point of interpolation of the x,y graph (please ensure
        that x is increasing order!), but only interpolate if the vector is
        longer than 12000 elements."""
        from scipy.interpolate import interp1d
        if len(x) >= 12000:
            x_new = nx.linspace(x[0], x[-1], 8000)
            interpolator = interp1d(x,y)
            y_new = interpolator(x_new)
            return x_new, y_new
        else:
            return x, y



    def plot_roc(self, fname, FPR, TPR, marker_FPR):
        """ROC curve (TPR vs FPR)"""
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        logging.debug("Plotting ROC curve to %s", fname.basename())
        # We get FPR and TPR decreasing, but interpolation is increasing-only
        FPR, TPR = self.interpolate(FPR[::-1], TPR[::-1]) 
        g.reset()
        g.title("ROC curve (TPR vs FPR)")
        g.ylabel("True Positive Rate (TPR)")
        g.xlabel("False Positive Rate (FPR)")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(FPR, TPR, title="TPR", with="lines"),
               Data([marker_FPR, marker_FPR], [0,0.99], title="threshold", with="lines"))


    def plot_precision(self, fname, TPR, PPV, marker_TPR):
        """Precision vs recall"""
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        logging.debug("Plotting Precision-Recall curve to %s", fname.basename())
        # We get TPR decreasing as threshold increases, reverse for interpolation
        TPR, PPV = self.interpolate(TPR[::-1], PPV[::-1])
        g.reset()
        g.title("Precision vs Recall")
        g.ylabel("Precision")
        g.xlabel("Recall")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(TPR, PPV, title="Precision", with="lines", smooth="csplines"),
               Data([marker_TPR, marker_TPR], [0,0.99], title="threshold", with="lines"))


    def plot_fmeasure(self, fname, pscores, TPR, PPV, FM, FMa, threshold):
        """Precision, Recall, F-Measure vs threshold"""
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        logging.debug("Plotting F-Measure curve to %s", fname.basename())
        npscores, TPR = self.interpolate(pscores, TPR)
        npscores, PPV = self.interpolate(pscores, PPV)
        npscores, FM  = self.interpolate(pscores, FM)
        npscores, FMa = self.interpolate(pscores, FMa)
        g.reset()
        g.title("Precision and Recall vs Threshold")
        g.ylabel("Precision, Recall, F-Measure, F-Measure Alpha")
        g.xlabel("Threshold Score")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(npscores, TPR, title="Recall",     with="lines"),
               Data(npscores, PPV, title="Precision",  with="lines"),
               Data(npscores, FM,  title="F1 Measure", with="lines"),
               Data(npscores, FMa, title="F Measure",  with="lines"),
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
        return min(150, max(10, bins))


    def plot_score_histogram(self, fname, pdata, ndata, threshold):
        """Histograms for pos and neg scores, with line to mark threshold""" 
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        logging.debug("Plotting article score histogram to %s", fname.basename())
        from itertools import chain
        py, px = nx.histogram(pdata, bins=self.bincount(pdata), normed=True)
        zy, zx = nx.histogram(ndata, bins=self.bincount(ndata), normed=True)
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


    def plot_feature_histogram(self, fname, scores):
        """Histogram for feature scores. Requires that there be no
        -1.#IND (undefined 0/0) or -1.#INF (zero-division like 1/0) values,
        as it messes up the histogram.  Only provide scores for
        features present in the inital corpus.  Features with DF=0
        (not present in a given corpus, although present somewhere in Medline) 
        far outnumber DF>0 features, making the histogram look like a spike.
        
        @param scores: List of feature weights."""
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        logging.debug("Plotting feature score histogram to %s", fname.basename())
        # Sorting modifies the vector, but the provided vector must be
        # left untouched to maintain the association with PubMed IDs.
        sscores = scores.copy() 
        sscores.sort()
        y, x = nx.histogram(scores, bins=self.bincount(sscores))
        y[y==0] = 1.0 # Remove infinities in the log axis
        g.reset()
        g.title("Feature Score Histogram")
        g.xlabel("Feature Score")
        g.ylabel("Histogram Mass")
        g("set logscale y")
        g("set terminal png")
        g("set output '%s'" % fname)
        g("set style fill solid 1.0")
        g.plot(Data(x, y, with="boxes"))



class DensityPlotter(Plotter):
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


    def plot_score_density(self, fname, pdata, ndata, threshold):
        """Probability density of pos and neg scores, with line to mark threshold
    
        @param pdata: Scores of positive documents
        @param ndata: Scores of negative documents
        @param threshold: Threshold score for counting a document positive
        """ 
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        from itertools import chain
        logging.debug("Plotting article score density to %s", fname.basename())
        px, py = self.gaussian_kernel_pdf(pdata)
        zx, zy = self.gaussian_kernel_pdf(ndata)
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


    def plot_feature_density(self, fname, scores):
        """Probability density function for feature scores"""
        if fname.exists() and not self.overwrite: return
        g = self.gnuplot
        logging.debug("Plotting feature score density to %s", fname.basename())
        x, y = self.gaussian_kernel_pdf(scores, npoints=1024)
        g.reset()
        g.title("Feature Score Density")
        g.xlabel("Feature Score")
        g.ylabel("Probability Density")
        g("set terminal png")
        g("set output '%s'" % fname)
        g.plot(Data(x, y, with="lines"))


