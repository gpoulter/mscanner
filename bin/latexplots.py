#!/usr/bin/env python

"""
Draws publication-quality plots for use in the paper

The figures are::
    Figure 1. Four score density plots
    Figure 2. ROC curve overlay
    Figure 3. PR curve overlay
    Figure 4. P,R,F1,Fa curve for AIDSBio to demo optimisation

"""

from __future__ import with_statement
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
from bin import retrievaltest


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


def read_featscores(indir, dataset):
    """Read feature scores for a dataset"""
    f = open(indir/dataset/mrc.report_term_scores, "r")
    f.readline()
    return array([float(s.split(",",1)[0]) for s in f])


def load_stats(indir, dataset, title, alpha=0.5):
    """Read statistics based on score data
    
    @param indir: Directory in which to find data sets
    @param dataset: Subdirectory name for the particular data set
    @param title: Name to put on the graphs (typically same as dataset)
    @param alpha: What alpha to use when recalcing performance
    """
    logging.info("Reading dataset %s", dataset)
    if not (indir/dataset).isdir():
        raise ValueError("Could not directory %s" % (indir/dataset)) 
    pscores = array([s[0] for s in scorefile.read_pmids(
        indir/dataset/mrc.report_positives, withscores=True)])
    nscores = array([s[0] for s in scorefile.read_pmids(
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


def calculate_overlap(px, py, nx, ny):
    """Calculate overlap between two bell curves (curve p must be
    to the right of curve n).  
    
    @deprecated: Overlap is not a meaningful performance statistic.
    
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


def plot_score_density(fname, statlist):
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
        # We don't plot overlapping area any more
        #area, iX, iY = plotting.calculate_overlap(px, py, nx, ny)
        #patch_overlap, = fill(iX, iY, facecolor='magenta', alpha=0.7, label=r"$\rm{Overlap}$")
        if idx == 0 or idx == 2:
            ylabel("Density")
        if idx == 2 or idx == 3:
            xlabel("Score")
        if idx == 0:
            legend(loc="upper left")
    custom_show(fname)


def plot_score_histogram(fname, pscores, nscores):
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


def plot_featscore_histogram(fname, fscores):
    """Plot histogram for individual feature scores"""
    logging.info("Plotting feature score histogram to %s", fname)
    ##title("Feature Score Histogram")
    xlabel("Feature Score")
    ylabel("Number of Features")
    fscores.sort()
    n, bins, patches = hist(fscores, bins=plotting.bincount(fscores))
    setp(patches, 'facecolor', 'r', 'linewidth', 0.0)
    custom_show(fname)


def plot_roc(fname, statlist):
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


def plot_precision(fname, statlist):
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


def plot_fmeasure(fname, s):
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

#### FUNCTIONS THAT USE THE ABOVE ####

def do_retrievaltest():
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
    pgx1_c = retrievaltest.compare_results_to_standard(pgx1, set(pg_test))
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


def do_publication():
    """Draws figures for the BMC paper: including densities, ROC curve, PR
    curve, and PRF curve. """
    indir = source_dir / "070622 CV10 100k a_i"
    aids = load_stats(indir, "aids-vs-100k", "AIDSBio")
    rad = load_stats(indir, "radiology-vs-100k", "Radiology")
    pg07 = load_stats(indir, "pg07-vs-100k", "PG07")
    ran10 = load_stats(indir, "random10k-vs-100k", "Random")
    all = (aids,rad,pg07,ran10)
    plot_roc("fig2_roc", all)
    plot_precision("fig3_pr", all)
    fmplot = validation.PerformanceStats(pg07.pscores, pg07.nscores, alpha=0.95)
    fmplot.title = pg07.title
    plot_fmeasure("fig6_prf", fmplot)
    plot_score_density("fig1_density", all)
    do_retrievaltest()


def do_testplots():
    """Tests the plot functions using some old smaller datasets"""
    global indir
    indir = source_dir / "Old Validation" / "070223 CV10 Daniel 2s"
    if not indir.isdir():
        raise ValueError("Cannot find %s" % indir)
    pg04 = load_stats(indir, "pg04-vs-30k", "PG04")
    pg07 = load_stats(indir, "pg07-vs-30k", "PG07")
    plot_score_density("test_density", (pg04,pg07,pg07,pg04))
    plot_roc("test_roc", (pg04,pg07))
    plot_precision("test_pr", (pg04,pg07))
    pg04alpha = validation.PerformanceStats(pg04.pscores, pg04.nscores, alpha=0.9)
    pg04alpha.title = pg04.title
    plot_fmeasure("test_prf", pg04alpha)


def do_subdirplots(subdirs):
    """Plots selected graphs for the datasets passed as parameters"""
    statlist = [load_stats(source_dir/path(d), d, d) for d in subdirs]
    fscores = [read_featscores(source_dir/path(d), d) for d in subdirs]
    #plot_score_density("cus_density", stats)
    #plot_roc("cus_roc", statlist)
    #plot_precision("cus_pr", statlist)
    for stats, fscores in zip(statlist, fscores):
        #plot_fmeasure("cus_%s_prf" % stats.title, stats)
        plot_score_histogram(
            "cus_%s_arthist" % stats.title, stats.pscores, stats.nscores)
        plot_featscore_histogram(
            "cus_%s_feathist" % stats.title, fscores)


if __name__ == "__main__":
    initLogger(logfile=False)
    if len(sys.argv) != 2:
        print "Please give python expression"
    else:
        eval(sys.argv[1])
