"""Cross-validation and performance statistic calculation"""

from __future__ import division

from itertools import chain, izip
import logging as log
import numpy as nx

from mscanner.FeatureScores import FeatureScores, FeatureCounts
from mscanner import utils


                                     
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


class Validator:
    """Cross-validated calculation of article scores.
    
    @group Constructor Parameters: featdb,featinfo,positives,negatives,nfolds,alpha,postfilter

    @ivar featdb: Mapping from doc id to list of feature ids
    
    @ivar featinfo: FeatureScores instance for stuff about features
    
    @ivar positives: Array of positive PMIDs for validation
    
    @ivar negatives: Array of negative PMIDs for validation
    
    @ivar nfolds: Number of validation folds (0 for leave-out-one)
    
    @group From validate: pscores,nscores
    
    @ivar pscores: Scores of positive articles after validation
    
    @ivar nscores: Scores of negative articles after validation
    """

    def __init__(self, featdb, featinfo,  positives,  negatives, nfolds = 10):
        """Constructor parameters set corresponding instance attributes."""
        pscores = None
        nscores = None
        utils.update(self, locals())


    def validate(self):
        """Carry out validation. 
        
        When L{nfolds}==0 we use leave-out-one instead of k-fold.
        
        @return: L{pscores}, L{nscores}"""
        if self.nfolds:
            return self.nfold_validate()
        else:
            return self.leaveout_validate()


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


    def nfold_validate(self, randomise=True):
        """Perform n-fold validation and return the raw performance measures
        
        @param randomise: Randomise validation splits (use False for debugging)
        
        @return: L{pscores}, L{nscores}
        """
        s = self
        pdocs = len(s.positives)
        ndocs = len(s.negatives)
        log.info("%d pos and %d neg articles", pdocs, ndocs)
        if randomise:
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
            log.debug("Fold %d: pstart = %d, psize = %s; nstart = %d, nsize = %d", 
                      fold, pstart, psize, nstart, nsize)
            # Get new feature scores
            s.featinfo.update_features(
                pos_counts = pcounts - FeatureCounts(
                    len(s.featinfo), s.featdb, 
                    s.positives[pstart:pstart+psize]), 
                neg_counts = ncounts - FeatureCounts(
                    len(s.featinfo), s.featdb, 
                    s.negatives[nstart:nstart+nsize]),
                pdocs = pdocs-psize, 
                ndocs = ndocs-nsize)
            # Calculate the article scores for the test fold
            termscores = s.featinfo.scores
            s.pscores[pstart:pstart+psize] = [
                nx.sum(termscores[s.featdb[d]]) for d in 
                s.positives[pstart:pstart+psize]]
            s.nscores[nstart:nstart+nsize] = [
                nx.sum(termscores[s.featdb[d]]) for d in 
                s.negatives[nstart:nstart+nsize]]
        return s.pscores, s.nscores


    def leaveout_validate(self):
        """Carries out leave-out-one validation, returning the resulting scores.
        
        @note: Feature scores by Bayesian pseudocount only - no other methods.
        
        @deprecated: 10-fold is standard, and leave-out-one is rather slow.
        
        @return: L{pscores}, L{nscores}
        """
        # Set up base feature scores
        pcounts = FeatureCounts(len(self.featinfo), self.featdb, self.positives)
        ncounts = FeatureCounts(len(self.featinfo), self.featdb, self.negatives)
        self.pscores = nx.zeros(len(self.positives), nx.float32)
        self.nscores = nx.zeros(len(self.negatives), nx.float32)
        pdocs = len(self.positives)
        ndocs = len(self.negatives)
        mask = self.featinfo.mask
        # Set up pseudocount
        if isinstance(self.featinfo.pseudocount, float):
            ps = nx.zeros(len(self.featinfo), nx.float32) + self.featinfo.pseudocount
        else:
            ps = self.featinfo.pseudocount
        marker = 0
        # Discount this article in feature score calculations
        def score_article(pmid, p_mod, n_mod):
            f = [fid for fid in self.featdb[doc] if not mask or not mask[fid]]
            return nx.sum(nx.log(
                ((pcounts[f]+p_mod+ps[f])/(pdocs+p_mod+2*ps[f]))/
                ((ncounts[f]+n_mod+ps[f])/(ndocs+n_mod+2*ps[f]))))
        # Get scores for positive articles
        for idx, doc in enumerate(self.positives):
            self.pscores[idx] = score_article(doc, -1, 0)
        # Get scores for negative articles
        for idx, doc in enumerate(self.negatives):
            self.nscores[idx] = score_article(doc, 0, -1)
        return self.pscores, self.nscores




