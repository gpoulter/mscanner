"""Provides a mapping between features and integer IDs"""

from __future__ import with_statement
import codecs
import numpy as nx


                                     
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


class FeatureMapping:
    """Persistent mapping between string features and feature IDs

    @note: Feature classes are "mesh", "qual", "issn", "word".  The class 
    identifies the source of the feature, and allows the same string
    to represent two different features - e.g. "human" as a "word" feature,
    and "human" as a "mesh" feature.
    
    @note: This is really a table with columns (id,type,name,count), and keys of id
    and (type,name).
    
    @group Passed via constructor: filename, ftype

    @ivar filename: Path to save the mapping to (None for memory only).
    
    @ivar ftype: Numpy type in feature vectors (uint16 or uint32).
    
    @ivar numdocs: Number of documents added while creating the mapping.
    
    @ivar features: Lookup list from feature ID to class and string
    (features[id] == (name,class)).
    
    @ivar feature_ids: Lookup from feature class and string to feature ID
    (feature_ids[class][name] == id)
    
    @ivar counts: Number of occurrences of each feature ID (counts[id] == #)
    in Medline.
    """

    def __init__(self, filename, ftype):
        """Initialise the database"""
        self.ftype = ftype
        self.filename = filename
        self.numdocs = 0
        self.features = []
        self.feature_ids = {}
        self.counts = []
        if filename is not None and self.filename.exists():
            self.load()


    def load(self):
        """Load feature mapping mapping from file as a tab-separated
        table for tuples (feature, type, count) with ID being the
        0-based index in the file.
        """
        self.features = []
        self.feature_ids = {}
        with codecs.open(self.filename, "rb", "utf-8") as f:
            self.numdocs = int(f.readline().strip())
            for fid, line in enumerate(f):
                feat, ftype, count = line.strip().split("\t")
                self.features.append((feat,ftype))
                self.counts.append(int(count))
                if ftype not in self.feature_ids:
                    self.feature_ids[ftype] = {}
                self.feature_ids[ftype][feat] = fid


    def dump(self):
        """Write the feature mapping to disk as a table of (feature, type,
        count) where feature ID+1 is the line number."""
        if self.filename is None:
            return
        _filename_new = self.filename + ".new"
        with codecs.open(_filename_new, "wb", "utf-8") as f:
            f.write("%s\n" % self.numdocs)
            for (feat, ftype), count in zip(self.features, self.counts):
                f.write("%s\t%s\t%d\n" % (feat, ftype, count))
        if self.filename.isfile():
            self.filename.remove()
        _filename_new.rename(self.filename)


    def __getitem__(self, key):
        """Get feature string, or feature ID depending on input.
        @param key: If feature ID return (feature, feature type). 
                    If (feature, feature type) returns feature ID
        """
        if isinstance(key, int):
            return self.features[key]
        elif isinstance(key, tuple) and len(key) == 2:
            return self.feature_ids[key[1]][key[0]]
        else:
            raise KeyError("Invalid key: %s" % str(key))


    def __len__(self):
        """Return number of distinct features"""
        return len(self.features)


    def class_mask(self, classes):
        """Boolean mask for features of particular classes.
        @param clases: Feature classes to exclude, like ["mesh","issn"].
        @return: Boolean array for features of those classes.
        """
        mask = nx.zeros(len(self.features), nx.bool)
        for fclass in classes:
            for fid in self.feature_ids[fclass].itervalues():
                mask[fid] = True
        return mask
    
    
    def get_vector(self, features):
        """Calculate the feature vector of an already-added article.
        
        @param features: Document features as a mapping from feature
        class to feature strings, like C{{'mesh':['A','B']}}.
     
        @return: Feature vector of corresponding feature IDs.
        """
        nfeats = sum(len(v) for v in features.itervalues())
        vector = nx.zeros(nfeats, self.ftype)
        idx = 0
        for fclass, featlist in features.iteritems():
            for feature in featlist:
                vector[idx] = self.feature_ids[fclass][feature]
                idx += 1
        # See FeatureStream - sorted vectors can be compressed.
        vector.sort() 
        return vector


    def add_article(self, features):
        """Add an article to the feature map.  Add 1 to number of 
        documents and to the count of each feature in the dictionary.
        Features not present in the map are created with a count of 1.
        
        @param features: Document features as a mapping from feature
        class to feature strings, like C{{'mesh':['A','B']}}."""
        self.numdocs += 1
        for fclass, featlist in features.iteritems():
            if fclass not in self.feature_ids:
                self.feature_ids[fclass] = {}
            fdict = self.feature_ids[fclass]
            for feat in featlist:
                if feat not in fdict:
                    featid = len(self.features)
                    self.features.append((feat,fclass))
                    self.counts.append(1)
                    fdict[feat] = featid
                else:
                    self.counts[fdict[feat]] += 1
