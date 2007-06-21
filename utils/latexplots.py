#!python

"""
Draws publication-quality plots for use in the paper

Figure 1. Four score density plots
Figure 2. ROC curve overlay
Figure 3. PR curve overlay
Figure 4. P,R,F1,Fa curve for AIDSBio to demo optimisation

"""

import logging
from path import path
from pylab import *

from mscanner.configuration import rc as mrc, initLogger
from mscanner.scorefile  import readPMIDs
from mscanner.validation import PerformanceStats
from mscanner.plotting import bincount, kernelPDF, calculateOverlap

initLogger()

interactive = False
npoints = 400
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
    f = file(indir/dataset/mrc.report_term_scores, "r")
    f.readline()
    return array([float(s.split(",",1)[0]) for s in f])

def getStats(indir, dataset, title, alpha=0.5):
    """Read statistics based on score data
    
    @param dataset: Subdirectory name for the data set
    @param title: Name to put on the graphs (typically same as dataset)
    @param alpha: What alpha to use when recalcing performance
    """
    logging.info("Reading dataset %s", dataset)
    pscores = array([s[0] for s in readPMIDs(indir/dataset/mrc.report_positives, withscores=True)])
    nscores = array([s[0] for s in readPMIDs(indir/dataset/mrc.report_negatives, withscores=True)])
    stats = PerformanceStats(pscores, nscores, alpha)
    stats.title = title
    return stats

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
        px, py = kernelPDF(s.pscores)
        nx, ny = kernelPDF(s.nscores)
        subplot(2,2,idx+1)
        title(s.title)
        line_pos, = plot(px, py, color='red', label=r"$\rm{Positive}$")
        line_neg, = plot(nx, ny, color='blue', label=r"$\rm{Negative}$")
        line_threshold = axvline(s.threshold, color='green', linewidth=1, 
                                 label=r"$\rm{Threshold}$")
        # No longer plotting the overlapping area
        #area, iX, iY = calculateOverlap(px, py, nx, ny)
        #patch_overlap, = fill(iX, iY, facecolor='magenta', alpha=0.7, 
        #                      label=r"$\rm{Overlap}$")
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
    title("Article Score Histograms")
    xlabel("Article Score")
    ylabel("Article Density")
    p_n, p_bins, p_patches = hist(pscores, bins=bincount(pscores), normed=True)
    n_n, n_bins, n_patches = hist(nscores, bins=bincount(nscores), normed=True)
    setp(p_patches, 'facecolor', 'r', 'alpha', 0.50, 'linewidth', 0.0)
    setp(n_patches, 'facecolor', 'b', 'alpha', 0.50, 'linewidth', 0.0)
    p_y = normpdf(p_bins, mean(pscores), std(pscores))
    n_y = normpdf(n_bins, mean(nscores), std(nscores))
    p_l = plot(p_bins, p_y, 'r--', label=r"$\rm{Positive}$")
    n_l = plot(n_bins, n_y, 'b--', label=r"$\rm{Negative}$")
    custom_show(fname, type="svg")
    
def plotFeatureScoreHistogram(fname, fscores):
    """Plot histogram for individual feature scores"""
    logging.info("Plotting feature score histogram to %s", fname)
    title("Feature Score Histogram")
    xlabel("Feature Score")
    ylabel("Number of Features")
    fscores.sort()
    n, bins, patches = hist(fscores, bins=bincount(fscores))
    setp(patches, 'facecolor', 'r', 'linewidth', 0.0)
    custom_show(fname)

def plotROC(fname, statlist):
    """Plots ROC curves overlayed"""
    logging.info("Plotting ROC grid to %s", fname)
    figure(figsize=(10,5))
    values = [smooth(s.TPR[::-1], s.FPR[::-1]) for s in statlist]
    # Plot complete ROC curve
    subplot(121)
    title(r"ROC curves")
    ylabel(r"True Positive Rate (Recall)")
    xlabel(r"False Positive Rate (1-Specificity)")
    lines = [plot(FPR, TPR)[0] for TPR, FPR in values]
    plot([s.tuned.FPR for s in statlist], [s.tuned.TPR for s in statlist], "kx", markeredgewidth=1)
    xlim(0.0, 1.0)
    ylim(0.0, 1.0)
    # Plot zoomed in ROC curve
    subplot(122)
    title(r"Magnified ROC")
    xlabel(r"False Positive Rate (1-Specificity)")
    lines = [plot(FPR, TPR)[0] for TPR, FPR in values]
    amount = 0.25
    legend(lines, [r"$\rm{"+s.title+r"}$" for s in statlist], "lower right")
    axis([0.0, amount, 1-amount, 1.0])
    custom_show(fname)

def plotPR(fname, statlist):
    """Plots PR curves overlayed"""
    logging.info("Plotting PR curve to %s", fname)
    title(r"Precision versus Recall")
    ylabel(r"$\rm{Precision}\ (\pi)$")
    xlabel(r"$\rm{Recall}\ (\rho)$")
    # Pairs of TPR and PPV vectors for plotting
    values = [smooth(s.TPR[::-1], s.PPV[::-1]) for s in statlist]
    lines = [plot(TPR, PPV)[0] for TPR, PPV in values]
    # Place X marks at threshold
    plot([s.tuned.TPR for s in statlist], [s.tuned.PPV for s in statlist], "kx", markeredgewidth=2)
    # Draw legends
    legend(lines, [r"$\rm{"+s.title+r"}$" for s in statlist], (0.1, 0.3))
    axis([0.0, 1.0, 0.0, 1.0])
    custom_show(fname)

def plotPRF(fname, s):
    """Plots a single Precision/Recall/F-Measure curve"""
    logging.info("Plotting PRF curve to %s", fname)
    title(s.title + " performance versus threshold")
    ylabel("Performance Measures")
    xlabel("Decision Threshold")
    x = linspace(s.uscores[0], s.uscores[-1], npoints)
    x, TPR = smooth(s.uscores, s.TPR, x)
    x, PPV = smooth(s.uscores, s.PPV, x)
    x, FM = smooth(s.uscores, s.FM, x)
    x, FMa = smooth(s.uscores, s.FMa, x)
    plot(x, TPR, label=r"$\rm{Recall}\ (\rho)$")
    plot(x, PPV, label=r"$\rm{Precision}\ (\pi)$")
    plot(x, FM, label=r"$F_1$")
    plot(x, FMa, label=r"$F_{\alpha}\ (\alpha=0.9)$")
    axvline(s.threshold, label=r"$\rm{Threshold}$")
    ylim(0,1)
    legend(loc="upper left")
    custom_show(fname)

def publication_plots():
    """Draws figures for the BMC paper: including densities, ROC curve, PR curve, and PRF curve.

    @note: Must still add the retrieval testing curve.
    """
    indir = path(r'C:\Documents and Settings\Graham\My Documents\data\mscanner\output\070530 CV10 vs 100k')
    aids = getStats(indir, "aids-vs-100k", "AIDSBio")
    rad = getStats(indir, "radiology-vs-100k", "Radiology")
    pg07 = getStats(indir, "pg07-vs-100k", "PG07")
    ran10 = getStats(indir, "random10k-vs-100k", "Random")
    all = (aids,rad,pg07,ran10)
    plotROC("fig2_roc", all)
    plotPR("fig3_pr", all)
    aidsalpha = PerformanceStats(aids.pscores, aids.nscores, alpha=0.95)
    aidsalpha.title = aids.title
    plotPRF("fig4_prf", aidsalpha)
    plotArticleScoreDensity("fig1_density", all)

def test_plots():
    global indir
    indir = path(r'C:\Documents and Settings\Graham\My Documents\data\mscanner\output\070206 LOO w s=0.01')
    pg04 = getStats(indir, "pg04-vs-30k", "PG04")
    pg07 = getStats(indir, "pg07-vs-30k", "PG07")
    plotArticleScoreDensity("test_density", (pg04,pg07,pg07,pg04))
    plotROC("test_roc", (pg04,pg07))
    plotPR("test_pr", (pg04,pg07))
    pg04alpha = PerformanceStats(pg04.pscores, pg04.nscores, alpha=0.9)
    pg04alpha.title = pg04.title
    plotPRF("test_prf", pg04alpha)

def custom_plots(subdirs):
    indir = path(r'C:\Documents and Settings\Graham\My Documents\data\mscanner\output\070206 LOO w s=0.01')
    statlist = [getStats(indir, d, d) for d in subdirs]
    fscores = [getFeatScores(indir, d) for d in subdirs]
    #plotArticleScoreDensity("cus_density", stats)
    #plotROC("cus_roc", statlist)
    #plotPR("cus_pr", statlist)
    for stats, fscores in zip(statlist, fscores):
        #plotPRF("cus_%s_prf" % stats.title, stats)
        plotArticleScoreHistogram("cus_%s_arthist" % stats.title, stats.pscores, stats.nscores)
        plotFeatureScoreHistogram("cus_%s_feathist" % stats.title, fscores)

if __name__ == "__main__":
    from getopt import getopt
    optlist, args = getopt(sys.argv[1:], "d:", ["data="])
    data = None
    for key, value in optlist:
        if key in ["-d","--data"]:
            data = value
    if data is None:
        print __doc__
        sys.exit(1)
    elif data== "publication":
        publication_plots()
    elif data == "test":
        test_plots()
    elif data == "custom":
        custom_plots(sys.argv[2:])
