"""Cross-validation and performance statistic calculation"""

from __future__ import division
import logging
import numpy as nx
import random

from mscanner import update
from mscanner.core.FeatureScores import FeatureScores


                                     
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


def make_partitions(nitems, nparts):
    """Partition a list of items into equal parts for cross validation.
    @param nitems: Number of items to partition.
    @param nparts: Number of partitions to form.
    @return: (starts, sizes) a list of start indeces and lengths of the partitions.
    """
    base, rem = divmod(nitems, nparts)
    sizes = base * nx.ones(nparts, nx.int32)
    sizes[:rem] += 1
    starts = nx.zeros(nparts, nx.int32)
    starts[1:] = nx.cumsum(sizes[:-1])
    return starts, sizes


def count_features(nfeats, featurevectors):
    """Count occurrenes of each feature in a set of articles.
    @param nfeats: Size of feature space.
    @param featurevectors: Iterable of feature vectors.
    @return: Array of length L{nfeats}, with the number of occurrences of each feature.
    """
    counts = nx.zeros(nfeats, nx.uint32)
    for featvec in featurevectors:
        counts[featvec] += 1
    return counts


def cross_validate(featinfo, positives, negatives, nfolds):
    """Perform cross validation. Remember to shuffle L{positives} and
    L{negatives} before passing them in!
    @param featinfo: L{FeatureScores} instance for getting feature scores.
    @param positives: List of feature vectors of relevant articles.
    @param negatives: List of feature vectors of irrelevant articles.
    @param nfolds: Number of validation folds
    @return: (pscores, nscores), vectors of document scores of the articles
    """
    pdocs = len(positives)
    ndocs = len(negatives)
    nfeats = len(featinfo.featmap)
    logging.debug("Cross-validating %d pos and %d neg items", pdocs, ndocs)
    pstarts, psizes = make_partitions(pdocs, nfolds)
    nstarts, nsizes = make_partitions(ndocs, nfolds)
    pscores = nx.zeros(pdocs, nx.float32)
    nscores = nx.zeros(ndocs, nx.float32)
    pcounts = count_features(nfeats, positives)
    ncounts = count_features(nfeats, negatives)
    for fold, (pstart,psize,nstart,nsize) in \
        enumerate(zip(pstarts,psizes,nstarts,nsizes)):
        logging.debug("Fold %d: pstart = %d, psize = %s; nstart = %d, nsize = %d", 
                  fold, pstart, psize, nstart, nsize)
        # Calculate feature scores for this fold
        featinfo.update(
            pos_counts = pcounts - count_features(nfeats, positives[pstart:pstart+psize]), 
            neg_counts = ncounts - count_features(nfeats, negatives[nstart:nstart+nsize]),
            pdocs = pdocs - psize, 
            ndocs = ndocs - nsize,
            prior = nx.log(pdocs/ndocs),
        )
        # Calculate the article scores for the test fold
        pscores[pstart:pstart+psize] =\
               featinfo.scores_of(positives[pstart:pstart+psize])
        nscores[nstart:nstart+nsize] =\
               featinfo.scores_of(negatives[nstart:nstart+nsize])
    return pscores, nscores

    
