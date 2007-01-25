from itertools import chain
import numpy
from Gnuplot import Data

def plotPDF(g, fname, pdata, ndata, threshold):
    """Plot PDFs for pos/neg scores, with line to mark threshold

    @param g: gnuplot object

    @param fname: Filename to plot to

    @param pdata: Scores of positive documents

    @param ndata: Scores of negative documents

    @param threshold: Threshold score for counting a document positive

    """ 
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
        
def plotHistograms(g, fname, pdata, ndata, threshold):
    """Plot histograms for pos/neg scores, with line to mark threshold
    (parameters as for plotPDF)
    """ 
    py, px = numpy.histogram(pdata)
    ny, nx = numpy.histogram(ndata)
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
