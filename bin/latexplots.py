#!/usr/bin/env python

"""
Draws publication-quality plots for use in the paper

The figures are::
    Figure 1. Four score density plots
    Figure 2. ROC curve overlay
    Figure 3. PR curve overlay
    Figure 4. P,R,F1,Fa curve for AIDSBio to demo optimisation

"""

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

import logging
from path import path

from pylab import *

from mscanner.configuration import rc as mrc, initLogger
from mscanner import plotting, scorefile, validation

interactive = False
npoints = 400
mscanner_dir = path(r"C:\Documents and Settings\Graham\My Documents\data\mscanner")
source_dir = mscanner_dir / "output"
outdir = path(r"C:\Documents and Settings\Graham\My Documents\temporary")

rc("figure", figsize=(8,6), dpi=100)
rc("figure.subplot", hspace=0.3)
rc("font", family="serif", serif="cmr10", monospace="cmtt10")
rc("legend", axespad=0.04, labelsep=0.015, pad=0.2)
rc("lines", linewidth=1.0)
rc("savefig", dpi=100)
rc("xtick.major", pad=6.0)
rc("ytick.major", pad=6.0)

def smooth(x, y, xn=npoints):
    """Resample a curve (x,y) using interpolation. xn is either a float with the
    new number of x value to evaluate the curve at, or an array of x values.
    """
    X = xn
    if not isinstance(X, arraytype):
        X = linspace(x[0], x[-1], xn)
    import scipy.interpolate
    interpolator = scipy.interpolate.interp1d(x, y, bounds_error=False)
    Y = interpolator(X)
    return X, Y

def getFeatScores(indir, dataset):
    """Read feature scores for a dataset"""
    f = open(indir/dataset/mrc.report_term_scores, "r")
    f.readline()
    return array([float(s.split(",",1)[0]) for s in f])

def getStats(indir, dataset, title, alpha=0.5):
    """Read statistics based on score data
    
    @param indir: Directory in which to find data sets
    @param dataset: Subdirectory name for the particular data set
    @param title: Name to put on the graphs (typically same as dataset)
    @param alpha: What alpha to use when recalcing performance
    """
    logging.info("Reading dataset %s", dataset)
    pscores = array([s[0] for s in scorefile.readPMIDs(
        indir/dataset/mrc.report_positives, withscores=True)])
    nscores = array([s[0] for s in scorefile.readPMIDs(
        indir/dataset/mrc.report_negatives, withscores=True)])
    stats = validation.PerformanceStats(pscores, nscores, alpha)
    stats.title = title
    return stats

def gplot(x,y,ls,label,pos=0.6,usemarker=False):
    """Wraps plot to add a single marker marker instead of lots"""
    if usemarker:
        i = int(pos*len(x))
        plot(x,y,ls[0:2])
        plot([x[i]],[y[i]],ls,label=label)
    else:
        plot(x,y,ls[0:2],label=label)

def custom_show(fname, doshow=interactive, type="eps"):
    """Either shows an interactive plot, or writes to EPS followed by
    convertion to PDF (since matplotlib's PDF backend is buggy)"""
    if interactive: 
        show()
    try:
        fullname = outdir/fname + "." + type
        savefig(fullname)
        if type == "eps":
            from subprocess import call
            call(["epstopdf", fullname], shell=True)
            import os
            os.remove(fullname)
    finally:
        close()

def plotArticleScoreDensity(fname, statlist):
    """Plots four score density plots in a grid

    @param statlist: A tuple of a four PerformanceStats objects
    """
    logging.info("Plotting score densities to %s", fname)
    for idx, s in enumerate(statlist):
        t = s.threshold
        px, py = plotting.kernelPDF(s.pscores)
        nx, ny = plotting.kernelPDF(s.nscores)
        subplot(2,2,idx+1)
        title(s.title)
        line_pos, = plot(px, py, color='red', label=r"$\rm{Positive}$")
        line_neg, = plot(nx, ny, color='blue', label=r"$\rm{Negative}$")
        line_threshold = axvline(
            s.threshold, color='green', linewidth=1, label=r"$\rm{Threshold}$")
        # No longer plotting the overlapping area
        #area, iX, iY = plotting.calculateOverlap(px, py, nx, ny)
        #patch_overlap, = fill(iX, iY, facecolor='magenta', alpha=0.7, label=r"$\rm{Overlap}$")
        if idx == 0 or idx == 2:
            ylabel("Density")
        if idx == 2 or idx == 3:
            xlabel("Score")
        if idx == 0:
            legend(loc="upper left")
    custom_show(fname)

def plotArticleScoreHistogram(fname, pscores, nscores):
    """Plot histograms for pos/neg scores, with line to mark threshold"""
    logging.info("Plotting article score histogram to %s", fname)
    ##title("Article Score Histograms")
    xlabel("Article Score")
    ylabel("Article Density")
    p_n, p_bins, p_patches = hist(pscores, bins=plotting.bincount(pscores), normed=True)
    n_n, n_bins, n_patches = hist(nscores, bins=plotting.bincount(nscores), normed=True)
    setp(p_patches, 'facecolor', 'r', 'alpha', 0.50, 'linewidth', 0.0)
    setp(n_patches, 'facecolor', 'b', 'alpha', 0.50, 'linewidth', 0.0)
    #p_y = normpdf(p_bins, mean(pscores), std(pscores))
    #n_y = normpdf(n_bins, mean(nscores), std(nscores))
    #p_l = plot(p_bins, p_y, 'r--', label=r"$\rm{Relevants}$")
    #n_l = plot(n_bins, n_y, 'b--', label=r"$\rm{Irrelevant}$")
    custom_show(fname)
    
def plotFeatureScoreHistogram(fname, fscores):
    """Plot histogram for individual feature scores"""
    logging.info("Plotting feature score histogram to %s", fname)
    ##title("Feature Score Histogram")
    xlabel("Feature Score")
    ylabel("Number of Features")
    fscores.sort()
    n, bins, patches = hist(fscores, bins=plotting.bincount(fscores))
    setp(patches, 'facecolor', 'r', 'linewidth', 0.0)
    custom_show(fname)

def plotROC(fname, statlist):
    """Plots ROC curves overlayed"""
    logging.info("Plotting ROC grid to %s", fname)
    figure(figsize=(10,5))
    formats = ["r-s", "b-D", "g-h", "c-"]
    values = [smooth(s.TPR[::-1], s.FPR[::-1]) for s in statlist]
    # Plot complete ROC curve
    subplot(121)
    ##title(r"ROC curves")
    ylabel(r"True Positive Rate (Recall)")
    xlabel(r"False Positive Rate (1-Specificity)")
    for (TPR, FPR), fmt  in zip(values, formats):
        plot(FPR, TPR, fmt[0:2])
    axis([0.0, 1.0, 0.0, 1.0])
    # Plot zoomed in ROC curve
    subplot(122)
    ##title(r"Magnified ROC")
    xlabel(r"False Positive Rate (1-Specificity)")
    for (TPR, FPR), fmt, s  in zip(values, formats, statlist):
        gplot(FPR, TPR, fmt, label=r"$\rm{"+s.title+r"}$", pos=0.96)
    amount = 0.25
    legend(loc="lower right")
    axis([0.0, amount, 1-amount, 1.0])
    custom_show(fname)

def plotPR(fname, statlist):
    """Plots PR curves overlayed"""
    logging.info("Plotting PR curve to %s", fname)
    ##title(r"Precision versus Recall")
    ylabel(r"$\rm{Precision}\ (\pi)$")
    xlabel(r"$\rm{Recall}\ (\rho)$")
    # Dotted line for break-even point
    plot([0.0, 1.0], [0.0, 1.0], "k:")
    # Pairs of TPR and PPV vectors for plotting
    formats = ["r-s", "b-D", "g-h", "c-o"]
    for s, fmt in zip(statlist, formats):
        TPR, PPV = smooth(s.TPR[::-1], s.PPV[::-1])
        gplot(TPR, PPV, fmt, label=r"$\rm{"+s.title+r"}$", pos=0.5)
    # Place X marks at threshold
    plot([s.tuned.TPR for s in statlist], 
         [s.tuned.PPV for s in statlist],
         "kx", markeredgewidth=2)
    # Draw legends
    legend(loc=(0.5, 0.15))
    axis([0.0, 1.0, 0.0, 1.0])
    custom_show(fname)

def plotPRF(fname, s):
    """Plots a single Precision/Recall/F-Measure curve"""
    logging.info("Plotting PRF curve to %s", fname)
    ##title(s.title + " performance versus threshold")
    ylabel("Performance Measures")
    xlabel("Decision Threshold")
    # We only plot starting from score of -50
    start = 0
    while s.uscores[start] < -40:
        start += 1
    x = linspace(s.uscores[start], s.uscores[-1], npoints)
    x, TPR = smooth(s.uscores[start:], s.TPR[start:], x)
    x, PPV = smooth(s.uscores[start:], s.PPV[start:], x)
    x, FM = smooth(s.uscores[start:], s.FM[start:], x)
    x, FMa = smooth(s.uscores[start:], s.FMa[start:], x)
    gplot(x, PPV, "b-D", label=r"$\rm{Precision}\ (\pi)$")
    gplot(x, FMa, "c-o", label=r"$F\ (\alpha=%s)$" % str(s.alpha))
    gplot(x, FM, "g-h", label=r"$F_1$")
    gplot(x, TPR, "r-s", label=r"$\rm{Recall}\ (\rho)$")
    axvline(s.threshold, c="k", label=r"$\rm{Threshold}$")
    ylim(0,1)
    legend(loc="upper right")
    custom_show(fname)

def RetrievalTest():
    """Plots retrieval test results for 20% of PharmGKB to see how MScanner
    and PubMed compare at retrieving the remaining 80%.
    """
    logging.info("Plotting Retrieval curve for PG07")
    pgdir = source_dir / "Retrieval" / "070626 Retrieval a_i 0.2 5k" / "pg07-retrieval"
    # Retrieval vs rank for MScanner
    mscanner_c = array([int(x) for x in (pgdir/"retrieval_stats.txt").lines()])
    # Gold standard citations for PG07
    pg_test = [int(x) for x in (pgdir/"retrieval_test.txt").lines()]
    # PubMed query output
    pubmed = mscanner_dir/"support"/"PubMed Pharmacogenetics"
    pgx1 = [int(x) for x in (pubmed/"pgx1.txt").lines()]
    # Retrieval vs rank for PubMed
    from mscanner.scoring import retrievalTest    
    pgx1_c = retrievalTest(pgx1, set(pg_test))
    # Plot the graph
    ax1 = subplot(111)
    ##title("PG07 Retrieval Comparison")
    ylabel("Fraction retrieved")
    xlabel("False Positives (logarithmic scale)")
    N = len(pg_test)
    r = 4301
    semilogx(range(1,r)-mscanner_c[1:r], mscanner_c[1:r]/N, "r", label="MScanner")
    semilogx(range(1,r)-pgx1_c[1:r], pgx1_c[1:r]/N, "b", label="PubMed")
    legend(loc="upper left")
    grid(True)
    gca().xaxis.grid(True, which='minor')
    ax2 = twinx()
    semilogx([1],[1])
    ax2.set_xlim(1,r)
    ax2.set_ylim(0,ax1.get_ylim()[1]*len(pg_test))
    ylabel("True Positives")
    custom_show("fig4_pg07retrieval")

def Publication():
    """Draws figures for the BMC paper: including densities, ROC curve, PR
    curve, and PRF curve.
    """
    indir = source_dir / "070621 CV10 100k a_i"
    aids = getStats(indir, "aids-vs-100k", "AIDSBio")
    rad = getStats(indir, "radiology-vs-100k", "Radiology")
    pg07 = getStats(indir, "pg07-vs-100k", "PG07")
    ran10 = getStats(indir, "random10k-vs-100k", "Random")
    all = (aids,rad,pg07,ran10)
    plotROC("fig2_roc", all)
    plotPR("fig3_pr", all)
    fmplot = validation.PerformanceStats(pg07.pscores, pg07.nscores, alpha=0.95)
    fmplot.title = pg07.title
    plotPRF("fig6_prf", fmplot)
    plotArticleScoreDensity("fig1_density", all)
    RetrievalTest()

def Testing():
    global indir
    indir = source_dir / "070206 LOO w s=0.01"
    pg04 = getStats(indir, "pg04-vs-30k", "PG04")
    pg07 = getStats(indir, "pg07-vs-30k", "PG07")
    plotArticleScoreDensity("test_density", (pg04,pg07,pg07,pg04))
    plotROC("test_roc", (pg04,pg07))
    plotPR("test_pr", (pg04,pg07))
    pg04alpha = validation.PerformanceStats(pg04.pscores, pg04.nscores, alpha=0.9)
    pg04alpha.title = pg04.title
    plotPRF("test_prf", pg04alpha)

def Custom(subdirs):
    statlist = [getStats(source_dir/d, d, d) for d in subdirs]
    fscores = [getFeatScores(source_dir/d, d) for d in subdirs]
    #plotArticleScoreDensity("cus_density", stats)
    #plotROC("cus_roc", statlist)
    #plotPR("cus_pr", statlist)
    for stats, fscores in zip(statlist, fscores):
        #plotPRF("cus_%s_prf" % stats.title, stats)
        plotArticleScoreHistogram(
            "cus_%s_arthist" % stats.title, stats.pscores, stats.nscores)
        plotFeatureScoreHistogram(
            "cus_%s_feathist" % stats.title, fscores)

if __name__ == "__main__":
    initLogger(logfile=False)
    if len(sys.argv) != 2:
        print "Please give python expression"
    else:
        eval(sys.argv[1])
