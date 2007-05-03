"""Implements the FeatureMapping type

FeatureMapping -- Mapping between string features and 16-bit feature IDs

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import codecs
import numpy as nx

class FeatureMapping:
    """Curates a database of string features, providing methods to map
    between a feature string and an integer feature ID.

    @note: The type of a feature could be 'mesh', 'qual', 'issn',
    'year', or 'author'.  A given feature string could have more than
    one type.

    @ivar numdocs: Number of documents used in creating the mapping

    @ivar features: List which effectively maps id:(feature,type)

    @ivar feature_ids: Dictionary to map type:{feature:id}

    @ivar counts: List mapping ID to number of occurrences
    """

    def __init__(self, featfile=None):
        """Initialise the database

        @param featfile: Path to text file with list of terms
        """
        self.featfile = featfile
        if self.featfile is not None:
            self.featfile_new = featfile+".new"
        self.numdocs = 0
        self.features = []
        self.feature_ids = {}
        self.counts = []
        if featfile is not None and self.featfile.exists():
            self.load()
        
    def load(self):
        """Load feature mapping mapping from file

        @note: file format is '%(feature)\t%(type)\t%(count)\n'
        """
        self.features = []
        self.feature_ids = {}
        f = codecs.open(self.featfile, "rb", "utf-8")
        self.numdocs = int(f.readline().strip())
        for fid, line in enumerate(f):
            feat, ftype, count = line.strip().split("\t")
            self.features.append((feat,ftype))
            self.counts.append(int(count))
            if ftype not in self.feature_ids:
                self.feature_ids[ftype] = {feat:fid}
            else:
                self.feature_ids[ftype][feat] = fid
        f.close()

    def dump(self):
        """Write the feature mapping to disk"""
        if self.featfile is None:
            return
        f = codecs.open(self.featfile_new, "wb", "utf-8")
        f.write("%s\n" % self.numdocs)
        for (feat, ftype), count in zip(self.features, self.counts):
            f.write(feat+"\t"+ftype+"\t"+str(count)+"\n")
        f.close()
        if self.featfile.isfile():
            self.featfile.remove()
        self.featfile_new.rename(self.featfile)

    def __getitem__(self, key):
        """Given a feature ID, returns (feature, feature type), given 
        (feature, feature type), returns feature ID"""
        if isinstance(key, int):
            return self.features[key]
        elif isinstance(key, tuple) and len(key) == 2:
            return self.feature_ids[key[1]][key[0]]
        else:
            raise KeyError("Invalid key: %s" % str(key))

    def __len__(self):
        """Return number of distinct features"""
        return len(self.features)

    def featureTypeMask(self, exclude_types):
        """Get a mask for excluded features

        @param exclude_types: Types of features to exclude

        @return: Boolean array for excluded features (but returns None if exclude_types is None)
        """
        if not exclude_types:
            return None
        exclude_feats = nx.zeros(len(self.features), nx.bool)
        for ftype in exclude_types:
            for fid in self.feature_ids[ftype].itervalues():
                exclude_feats[fid] = True
        return exclude_feats

    def addArticle(self, **kwargs):
        """Add an article given a list of feature strings of with keyword
        arguments for feature types.

        @note: Dynamically creates new features IDs and feature types as necessary.

        @note: Takes keyword arguments for feature types, with values being a
        list of feature strings for that type. e.g. mesh=["Term A","Term B"]

        @return: Numpy array of uint16 feature IDs
        """
        result = []
        self.numdocs += 1
        for ftype, fstrings in kwargs.iteritems():
            if ftype not in self.feature_ids:
                self.feature_ids[ftype] = {}
            fdict = self.feature_ids[ftype]
            for feat in fstrings:
                if feat not in fdict:
                    featid = len(self.features)
                    self.features.append((feat,ftype))
                    self.counts.append(1)
                    fdict[feat] = featid
                else:
                    self.counts[fdict[feat]] += 1
                result.append(fdict[feat])
        return nx.array(result, nx.uint16)
