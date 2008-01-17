"""Cross-validation and performance statistic calculation"""

from __future__ import division
import logging
import numpy as nx

from mscanner import update
from mscanner.core.FeatureScores import FeatureScores, FeatureCounts


                                     
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


class CrossValidator:
    """Cross-validated calculation of article scores.
    
    @group Constructor Parameters: featdb,featinfo,positives,negatives,nfolds,alpha,postfilter

    @ivar featdb: Mapping from doc id to list of feature ids
    
    @ivar featinfo: L{FeatureScores} instance to handle training
    
    @ivar positives: Array of positive PMIDs for validation
    
    @ivar negatives: Array of negative PMIDs for validation
    
    @ivar nfolds: Number of validation folds
    
    @ivar randseed: Used with numpy.random.seed to get the same shuffle of the
    data before partitioning (to see how different parameters affect
    performance, without stochastic variation from the shuffle making it harder
    to tell how big the difference is).
    
    @group Set by validate: pscores,nscores
    
    @ivar pscores: Scores of positive articles after validation
    
    @ivar nscores: Scores of negative articles after validation
    """

    def __init__(self, featdb, featinfo, positives, negatives, nfolds, randseed=None):
        update(self, locals())


    @staticmethod
    def make_partitions(nitems, nparts):
        """Calculate partitions of input data for cross validation
        
        @param nitems: Number of items to partition
        @param nparts: Number of partitions
        @return: List of start indeces, and list of lengths for partitions
        """
        base, rem = divmod(nitems, nparts)
        sizes = base * nx.ones(nparts, nx.int32)
        sizes[:rem] += 1
        starts = nx.zeros(nparts, nx.int32)
        starts[1:] = nx.cumsum(sizes[:-1])
        return starts, sizes


    def validate(self, _randomise=True):
        """Perform n-fold validation and return the raw performance measures
        
        @param _randomise: If True, shuffle the corpora before partitioning
        into folds. Only use False for unit testing.
       
        @return: L{pscores}, L{nscores}
        """
        s = self
        pdocs = len(s.positives)
        ndocs = len(s.negatives)
        logging.debug("Cross-validating %d pos and %d neg items", pdocs, ndocs)
        if _randomise:
            logging.debug("Shuffling corpora with seed %s", str(s.randseed))
            nx.random.seed(s.randseed)
            nx.random.shuffle(s.positives)
            nx.random.shuffle(s.negatives)
        s.pstarts, s.psizes = s.make_partitions(pdocs, s.nfolds)
        s.nstarts, s.nsizes = s.make_partitions(ndocs, s.nfolds)
        s.pscores = nx.zeros(pdocs, nx.float32)
        s.nscores = nx.zeros(ndocs, nx.float32)
        pcounts = FeatureCounts(len(s.featinfo), s.featdb, s.positives)
        ncounts = FeatureCounts(len(s.featinfo), s.featdb, s.negatives)
        for fold, (pstart,psize,nstart,nsize) in \
            enumerate(zip(s.pstarts,s.psizes,s.nstarts,s.nsizes)):
            logging.debug("Fold %d: pstart = %d, psize = %s; nstart = %d, nsize = %d", 
                      fold, pstart, psize, nstart, nsize)
            # Get new feature scores
            s.featinfo.update(
                pos_counts = pcounts - FeatureCounts(
                    len(s.featinfo), s.featdb, 
                    s.positives[pstart:pstart+psize]), 
                neg_counts = ncounts - FeatureCounts(
                    len(s.featinfo), s.featdb, 
                    s.negatives[nstart:nstart+nsize]),
                pdocs = pdocs-psize, 
                ndocs = ndocs-nsize,
                prior = nx.log(pdocs/ndocs),
            )
            # Calculate the article scores for the test fold
            s.pscores[pstart:pstart+psize] = s.featinfo.scores_of(
                s.featdb, s.positives[pstart:pstart+psize])
            s.nscores[nstart:nstart+nsize] = s.featinfo.scores_of(
                s.featdb, s.negatives[nstart:nstart+nsize])
        return s.pscores, s.nscores
