from itertools import chain
import numpy as n
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from Gnuplot import Data

def calculateOverlap(px, py, nx, ny):
    """Given bell curves for positives and negatives, calculate the area of
    overlap, assuming that p is located to the right of n. 
    
    Procedure is first to interpolate both curves onto a common X-axis ranging
    from min(nx) to max(px), then find highest where pY and nY intersect,
    then create a curve using pY up to intersection, and nY thereafter. Then
    calculate the area underneath using trapezoidal rule."""
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

def plotPDF(g, fname, pdata, ndata, threshold):
    """Plot PDFs for pos/neg scores, with line to mark threshold

    @param g: gnuplot object
    @param fname: Filename to plot to
    @param pdata: Scores of positive documents
    @param ndata: Scores of negative documents
    @param threshold: Threshold score for counting a document positive
    """ 
    px, py = kernelPDF(pdata)
    nx, ny = kernelPDF(ndata)
    overlap = calculateOverlap(px, py, nx, ny)
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
    return overlap
        
def plotHistograms(g, fname, pdata, ndata, threshold):
    """Plot histograms for pos/neg scores, with line to mark threshold
    (parameters as for plotPDF)
    """ 
    py, px = n.histogram(pdata)
    ny, nx = n.histogram(ndata)
    g.ylabel("Number of articles")
    g.xlabel("Score of article")
    g.title("Score Histograms")
    g("set terminal png")
    g("set output '%s'" % fname)
    threshold_height = max(chain(py, ny))
    g.plot(Data([threshold, threshold], [0, threshold_height], title="threshold", with="lines"),
           Data(px, py, title="Positives", with="histeps"),
           Data(nx, ny, title="Negatives", with="histeps"))

def plotROC(g, roc, FPR, TPR, marker_FPR):
    """ROC curve (TPR vs FPR)

    @param g: gnuplot object

    @param roc: Path to output file for ROC curve
    """
    g.reset()
    g.ylabel("True Positive Rate (TPR)")
    g.xlabel("False Positive Rate (FPR)")
    g.title("ROC curve (TPR vs FPR)")
    g("set terminal png")
    g("set output '%s'" % roc)
    g.plot(Data(FPR, TPR, title="TPR", with="lines"),
           Data([marker_FPR, marker_FPR], [0,0.99], title="threshold", with="lines"))

def plotPrecisionRecall(g, p_vs_r, TPR, PPV, marker_TPR):
    """Precision vs recall graph

    @param g: gnuplot object

    @param p_vs_r: Path to output file for precision-recall curve
    """
    g.reset()
    g.ylabel("Precision")
    g.xlabel("Recall")
    g.title("Precision vs Recall")
    g("set terminal png")
    g("set output '%s'" % p_vs_r)
    g.plot(Data(TPR, PPV, title="Precision", with="lines", smooth="csplines"),
           Data([marker_TPR, marker_TPR], [0,0.99], title="threshold", with="lines"))

def plotPrecisionRecallFmeasure(g, pr_vs_score, pscores, TPR, PPV, FM, FMa, threshold):
    """Precision, Recall, F-Measure vs threshold graph

    @param g: gnuplot object

    @parav pr_vs_score: Path to output file for precision,recall vs threshold
    """
    g.reset()
    g.ylabel("Precision, Recall, F-Measure, F-Measure Alpha")
    g.xlabel("Threshold Score")
    g.title("Precision and Recall vs Threshold")
    g("set terminal png")
    g("set output '%s'" % pr_vs_score)
    g.plot(Data(pscores, TPR, title="Recall", with="lines"),
           Data(pscores, PPV, title="Precision", with="lines"),
           Data(pscores, FM, title="F1 Measure", with="lines"),
           Data(pscores, FMa, title="F Measure", with="lines"),
           Data([threshold, threshold], [0,0.99], title="threshold", with="lines"))
