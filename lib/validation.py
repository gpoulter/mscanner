#!env python

"""Cross-validation and performance statistics module

@author: Graham Poulter
                                   

Validator -- Perform cross-validation and output performance statistics

"""

from __future__ import division
import sys
import unittest
import logging as log
from random import randint, seed
from heapq import heapreplace
from itertools import chain
from path import path

from article import TermCounts
from scoring import getTermScores, scoreDocument, writeTermScoresCSV, writeTermScoresHTML
import templates

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
        recall,
        nfold,
        pseudocount=0.1,
        daniel=False,
        genedrug_articles=None ):
        """Initialise validator
        @param meshdb: Mapping termid:term (FeatureMapping), for writing term names in report
        @param featdb: docid:[termid] mapping
        @param posids: Set of positive docids
        @param negids: Set of negative docids
        @param recall: Calibrate threshold to achieve this recall on training data
        @param nfold: Number of folds in cross-validation.
        @param pseudocount, daniel: Passed to Scoring.getTermScores
        @param genedrug_articles: Set of docids having gene-drug co-occurrences
        """
        self.meshdb = meshdb
        self.featdb = featdb
        self.pos = posids
        self.neg = negids
        self.recall = recall
        self.nfold = nfold
        self.pseudocount = pseudocount
        self.daniel = daniel
        self.genedrug_articles = genedrug_articles

    def doFold( self, ptrain, ptest, ntrain, ntest ):
        """Carry out one inner iteration of the validation loop
        @param ptrain,ptest,ntrain,ntest: Positive and negative training / testing data sets
        @return: Return TP,TN,FP,FN, threshold, threshold, positive scores, negative scores
        """
        featdb = self.featdb
        pfreqs = TermCounts( featdb[d] for d in ptrain )
        nfreqs = TermCounts( featdb[d] for d in ntrain )
        termscores = getTermScores( pfreqs, nfreqs, self.pseudocount, self.daniel )
        if self.daniel:
            threshold = 10 # Hardcoded for Daniel's scoring method
        else:
            ptrain_scores = [ (d,scoreDocument(featdb[d],termscores)) for d in ptrain ]
            threshold = self.getThreshold( [s for d,s in ptrain_scores], self.recall )
        ptest_scores = [ (d,scoreDocument(featdb[d],termscores)) for d in ptest ]
        ntest_scores = [ (d,scoreDocument(featdb[d],termscores)) for d in ntest ]
        # Evaluate positive test articles to get TPs and FNs
        TP = 0
        for docid, score in ptest_scores:
            if self.articleIsPositive( docid, score, threshold ):
                TP += 1
        FN = len(ptest) - TP
        # Evaluate negative test articles to get FPs and TNs
        FP = 0
        for docid, score in ntest_scores:
            if self.articleIsPositive( docid, score, threshold ):
                FP += 1
        TN = len(ntest) - FP
        # Return overall performance statistics
        return TP, TN, FP, FN, threshold, [s for d,s in ptest_scores], [s for d,s in ntest_scores]

    def articleIsPositive( self, docid, score, threshold ):
        """Classifies an article as positive or negative based on score threshold"""
        if self.genedrug_articles is None or docid in self.genedrug_articles:
            return score >= threshold
        else:
            return False

    @staticmethod
    def getThreshold( scores, recall ):
        """Return score threshold for given recall
        @param scores: Log-likelihood scores of articles
        @param recall: Proportion of articles to score higher than return value
        @return: Threshold that achieves given recall
        """
        if recall < 0 or recall > 1:
            raise ValueError("recall %s is not between 0.0 and 1.0" % str(recall))
        scores.sort()
        cutindex = int( (1.0-recall) * len( scores ) )
        return scores[cutindex]

    @staticmethod
    def partition( items, nfold ):
        """Partitions indeces into data items for cross-validation

        @param items: Set of objects to be partitioned
        @param nfold: Number of partitions to use
        @return: List of sets, each set being length len(items)/nfold
        """
        tests = [ set() for n in range(nfold) ]
        ct = 0 # current test
        available = list(items) # items to assign
        N = len(available) # number of items
        i = 0 # current iteration
        size = N # number of items left
        while i < N:
            r = randint(0,size-1) # choose a random item
            tests[ct].add(available[r]) # add item to current set
            available[r] = available[size-1] # over write with last item
            i += 1 # next iteration
            ct = (ct+1) % len(tests) # next test
            size = size - 1 # virtual deletion
        return tests

    def validate(self):
        """Perform n-fold validation and return the raw performance measures"""
        positives = self.pos
        negatives = self.neg
        log.info( "%d pos and %d neg articles", len(positives), len(negatives) )
        TP,TN,FP,FN = (0,0,0,0)
        ptests = self.partition(positives, self.nfold)
        ntests = self.partition(negatives, self.nfold)
        pscores, nscores, threshold = [],[],0
        for fold, (ptest, ntest) in enumerate(zip(ptests, ntests)):
            ptrain = positives - ptest
            ntrain = negatives - ntest
            _TP,_TN,_FP,_FN,_threshold,_pscores,_nscores = self.doFold(ptrain, ptest, ntrain, ntest)
            pscores.extend( _pscores )
            nscores.extend( _nscores )
            log.info( "Fold %d: TP=%d, TN=%d, FP=%d, FN=%d", fold, _TP, _TN, _FP, _FN )
            TP,TN,FP,FN = ( TP+_TP, TN+_TN, FP+_FP, FN+_FN )
            threshold += _threshold
        recall = 0
        precision = 0
        fmeasure = 0
        if TP+FN != 0: recall = TP/(TP+FN)
        if TP+FP != 0: precision = TP/(TP+FP)
        if recall > 0 and precision > 0: fmeasure = 2*recall*precision/(recall+precision)
        threshold = threshold / self.nfold
        log.info( "Average: threshold=%g, recall=%g, precision=%g, fmeasure=%g",
                  threshold, recall, precision, fmeasure )
        return TP, TN, FP, FN, threshold, pscores, nscores

    @staticmethod
    def plotHistogramsPYLAB(name,pdata,ndata,threshold,recall,precision):
        """Plot histograms of the positive and negative scores, with a
        vertical line marking the classifier threshold."""
        import pylab as p
        binCount=15
        p_n,p_bins,p_patches=p.hist(pdata, bins=binCount, normed=1, alpha=1.0, facecolor='r')
        y = p.normpdf(p_bins, p.mean(pdata), p.std(pdata))
        l = p.plot(p_bins, y, 'r--', linewidth=1)
        n_n,n_bins,n_patches=p.hist(ndata, bins=binCount, normed=1, alpha=0.5, facecolor='b')
        y = p.normpdf(n_bins, p.mean(ndata), p.std(ndata))
        l = p.plot(n_bins, y, 'b--', linewidth=1)
        l = p.plot([threshold, threshold], [0, max(chain(p_n,n_n))], 'k', linewidth=3)
        lcolors=(p_patches[0], n_patches[0])
        p.xlabel('Article score')
        p.ylabel('Probability density')
        p.legend(lcolors,('Positive Set','Negative Set'))
        p.title('Score Densities (r=%d,p=%d)' % (int(recall*100),int(precision*100)))
        p.savefig(name)
        p.close()

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
    def plotHistograms(name, pdata, ndata, threshold):
        """Plot histograms of the positive and negative scores, with a
        vertical line marking the classifier threshold."""
        from Gnuplot import Gnuplot, Data
        g = Gnuplot(debug=1)
        g.xlabel('Article score')
        g.ylabel('Probability density')
        g.title('Score Densities')
        g('set terminal png')
        g('set output "'+name+'"')
        (pcen,pfreq) = Validator.makeHistogram(pdata)
        (ncen,nfreq) = Validator.makeHistogram(ndata)
        threshold_height = max( chain( pfreq, nfreq ) )
        g.plot(Data([threshold,threshold],[0,threshold_height],title="threshold",with='lines'),
               Data(pcen,pfreq,title='Positives',with='histeps'),
               Data(ncen,nfreq,title='Negatives',with='histeps'))
        
    @staticmethod
    def plotCurves(roc, p_vs_r, pr_vs_score, pscores, nscores):
        """Plot precision-recall, and precision,recall vs threshold"""

        # Calculate TPR (recall), FPR, and precision
        pscores.sort()
        nscores.sort()
        P = len(pscores)
        N = len(nscores)
        TPR = [ 0 for xi in xrange(P) ] # recall
        FPR = [ 0 for xi in xrange(P) ] # 1-specificity
        precision = [ 0 for xi in xrange(P) ] # precision
        TN = 1
        for xi in xrange(P):
            threshold = pscores[xi]
            while (nscores[TN-1] < threshold) and (TN < N):
                TN += 1
            FN = xi+1      # xi+1 is # of positives we called negative
            TP = P - FN    # TP+FN = P
            FP = N - TN    # TN+FP = N
            TPR[xi] = TP/P # TPR = TP/P
            FPR[xi] = FP/N # FPR = FP/N = 1 - TN/N = 1 - specificity
            if TP+FP > 0:
                precision[xi] = TP/(TP+FP) # PPV = TP/(TP+FP) = precision
            #print "thresh = %g, TPR = %d/%d = %.1e, FPR = %d/%d = %.1e" % (threshold, TP, P, TP/P, FP, N, FP/N)
        recall = TPR
            
        from Gnuplot import Gnuplot, Data
        g = Gnuplot(debug=1)

        # ROC (TPR vs FPR)
        g.xlabel("False Positive Rate (FPR)")
        g.ylabel("True Positive Rate (TPR)")
        g.title("ROC curve (TPR vs FPR)")
        g('set terminal png')
        g('set output "%s"' % roc)
        g.plot( Data( FPR, TPR, title="TPR", with="lines" ) )

        # Precision vs recall
        g.xlabel("Recall")
        g.ylabel("Precision")
        g.title("Precision vs Recall")
        g('set terminal png')
        g('set output "%s"' % p_vs_r)
        g.plot( Data( recall, precision, title="Precision", with="lines" ) )

        # Precision and recall vs threshold
        g.reset()
        g.xlabel('Threshold Score')
        g.ylabel('Precision and Recall')
        g.title('Precision and Recall vs Threshold')
        g('set terminal png')
        g('set output "%s"' % pr_vs_score)
        g.plot( Data( pscores, recall, title="Recall", with="lines" ),
                Data( pscores, precision, title="Precision", with="lines" ) )

    def report(self, results, prefix, stylesheet):
        """Write a full validation report

        @param resuls: (TP,TN,FP,FN,threshold,pscores,nscores) from self.validate()
        @param prefix: Prefix for output files of validation report
        @param stylesheet: Path to CSS stylesheet
        """
        # Set up output files
        terms_csv = prefix/"termscores.csv"
        terms_html = prefix/"termscores.html"
        hist_img = prefix/"histogram.png"
        roc_img = prefix/"roc.png"
        p_vs_r_img = prefix/"prcurve.png"
        pr_vs_score_img = prefix/"prscore.png"
        mainfile = prefix/"index.html"
        # Calculate performance measures
        TP, TN, FP, FN, threshold, pscores, nscores = results
        P = TP+FN
        N = TN+FP
        A = TP+FN+TN+FP # = P + N = T + F
        T = TP+TN
        F = FP+FN
        TPR, FNR, TNR, FPR, PPV, NPV = 0, 0, 0, 0, 0, 0
        accuracy, prevalence, fmeasure = 0, 0, 0
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
        # Term scores
        pfreqs = TermCounts(self.featdb[d] for d in self.pos)
        nfreqs = TermCounts(self.featdb[d] for d in self.neg)
        termscores = getTermScores( pfreqs, nfreqs, self.pseudocount, self.daniel)
        writeTermScoresCSV(file(terms_csv,"w"), self.meshdb, termscores, pfreqs, nfreqs)
        writeTermScoresHTML(file(terms_html,"w"), self.meshdb, termscores, pfreqs, nfreqs, self.pseudocount)
        # Output performance statistics
        self.plotHistograms(hist_img, pscores, nscores, threshold)
        self.plotCurves(roc_img, p_vs_r_img, pr_vs_score_img, pscores, nscores)
        templates.validation.run(dict(
            TP=TP, TN=TN, FP=FP, FN=FN, P=P, N=N, A=A, T=T, F=F,
            TPR=TPR, FNR=FNR, TNR=TNR, FPR=FPR, PPV=PPV, NPV=NPV,
            accuracy = accuracy,
            prevalence = prevalence,
            recall = recall,
            trainrecall = self.recall,
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
        
class _ValidatorTest(unittest.TestCase):
    def test(self):
        val = Validator(
            meshdb = { 1:"A", 2:"B", 3:"C", 4:"D", 5:"E", 6:"F" },
            featdb = { 1:[1,2,3], 2:[2,3], 3:[1,3], 4:[4,6], 5:[4], 6:[3,4,5], 7:[4,5] },
            posids = set([ 1, 2, 3 ]),
            negids = set([ 4, 5, 6, 7]),
            recall = 0.9,
            nfold = 2,
            pseudocount = 0.1,
            daniel = False,
            genedrug_articles = None,
            )
        seed(0)
        self.assertEqual(
            val.partition(set([1,2,3,4,5,6,7,8,9,10]),3) ,
            [set([9,2,10,6]),set([8,1,7]),set([3,4,5])] )
        results = val.validate()
        val.report(results, path("/tmp/valid_"), path("templates/style.css"))

if __name__ == "__main__":
    unittest.main()
