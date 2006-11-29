"""Configuration for MedScanner

@author: Graham Poulter
                                   

Import this module to set up library paths and configure locations to
files used by the scanner, as well as configuration options.

"""

import logging
import sys
import os
from new import module
from path import path

#### Logging configuration

logging.basicConfig(
    level    = logging.DEBUG,
    datefmt  = "%H:%M:%S",
    format   = "%(asctime)-9s %(levelname)-8s %(message)s",
)

#### TOP-LEVEL CONFIGURATION (DO NOT MODIFY FROM OUTSIDE)

## Path to source directory
src = path(sys.argv[0]).dirname().abspath().parent
## Path to top 
top = src.parent
## Path to lib
lib = src / "lib"
sys.path.insert( 0, lib )
## Path to templates
templates = lib / "templates"
## Path to base for working data
working = top / "data"
#working = top / "testing"
## Path to output files
output = working / "output"
## Path to cache directory
cache = working / "cache"
## Path to web output directory
weboutput = src / "www" / "htdocs" / "output"
## Path to corpora
corpora = working / "corpora"

#### DATABASE CONFIGURATION

## Path to stylesheet for query/validation reports
stylesheet = templates / "style.css"
## Path to file with progress statistics
statfile = cache / "mscanner.pid"
## Path to BSD DB environment home directory
db_home = cache / "db_home"
## Path to DB of article objects
articledb = cache / "articles.db"
## Path to list of article IDs
articlelist = cache / "articles.txt"
## Path to DB of term features for each article
featuredb = cache / "features.db"
## Logical database to use withing featuredb
featureset = "meshterms" # "allfeats", "meshterms"
## Path to feature<->ID mapping
featuremap = cache / "featuremap.txt"
## Path to term counts pickle
termcounts = cache / "termcounts.pickle"
## Path to list files already processed
processed = cache / "processed.txt"
## Whether to use transactions while updating
use_transactions = False

#### GENEDRUG CONFIGURATION

## Path to GAPScore results in
gapscore = cache / "gapscore.db"
## Path to gene-drug co-occurrence cache
genedrug = cache / "genedrug.db"
## Path to pickled drug table (generated by drugtable.py)
drugtable = working / "drug-table" / "drugtable.pickle"

#### UPDATE CONFIGURATION (for update.py)

## Path to compressed MEDLINE 
medline = working / "medline"
## Number of files to process between saving
save_delay = 3
## Pickle with MeSH terms to exlude (generated by meshexcludes.py)
mesh_excludes = working / "mesh-excludes" / "meshexcludes.pickle"
## Pickle with MeSH synonyms (generated by meshexcludes.py)
mesh_synonyms = working / "mesh-excludes" / "meshsynonyms.pickle"

#### QUERY CONFIGURATION (for query.py)

## Positive set
posfile = None
## Integer for maximum number of results (may be fewer due to threshold)
limit = 10000
## Float for minimum score threshold
threshold = 10
## Per-term pseudocount to use
pseudocount = 0.01
## Prefix for result report files
query_report = output / "results"
## Path to database output
outputdb = None
## Types of features to exclude
exclude_feats = None # ["issn"]

#### VALIDATOR CONFIGURATION (for validate.py)

## Positive set
posfile = None
## Negative set
negfile = None
## Validation folds to use
nfolds = 5
## Whether to do a gene-drug co-occurrence filter
dogenedrug = False
## Whether to use Daniel's 10^-8 pseudocounts
dodaniel = False
## Prefix to use for report file
valid_report = output / "validation"


#### Configure which data sets to ues

#dataset = "aids-vs-500k"
#dataset = "aidsbig-vs-500k"
#dataset = "cur-vs-500k"
#dataset = "cur-vs-go4"
#dataset = "cur-vs-med"
#dataset = "daniel-vs-500k"
#dataset = "full-vs-100k"
#dataset = "full-vs-500k"
#dataset = "full-vs-go4"
#dataset = "gdtest"
#dataset = "old-vs-30k"
#dataset = "old-vs-500k"
dataset = "old-vs-go4"

query_report = output / (dataset+"-query")
valid_report = output / (dataset+"-valid")

if dataset == "aids-vs-500k":
    posfile = corpora / "ncbi-aids-bioethics-Oct06-1k.txt"
    negfile = corpora / "medline-500k.txt"
if dataset == "aidsbig-vs-500k":
    posfile = corpora / "ncbi-aids-bioethics-Oct06.txt"
    negfile = corpora / "medline-500k.txt"
if dataset == "cur-vs-500k":
    posfile = corpora / "pharmgkb-Oct06.txt"
    negfile = corpora / "medline-500k.txt"
if dataset == "cur-vs-go4":
    posfile = corpora / "pharmgkb-Oct06.txt"
    negfile = corpora / "geneontology-2004.txt"
if dataset == "daniel-vs-500k":
    posfile = corpora / "daniel-small.txt"
    negfile = corpora / "medline-500k.txt"
if dataset == "cur-vs-med":
    posfile = corpora / "pharmgkb-Oct06.txt"
    negfile = articlelist
if dataset == "full-vs-500k":
    posfile = corpora / "pharmgkb-full-Jan06.txt"
    negfile = corpora / "medline-500k.txt"
if dataset == "full-vs-go4":
    posfile = corpora / "pharmgkb-full-Jan06.txt"
    negfile = corpora / "geneontology-2004.txt"
if dataset == "gdtest":
    posfile = corpora / "genedrug-small.txt"
    negfile = corpora / "geneontology-2004.txt"
if dataset == "old-vs-30k":
    posfile = corpora / "pharmgkb-2004.txt"
    negfile = corpora / "medline-30k.txt"
if dataset == "old-vs-500k":
    posfile = corpora / "pharmgkb-2004.txt"
    negfile = corpora / "medline-500k.txt"
if dataset == "old-vs-go4":
    posfile = corpora / "pharmgkb-2004.txt"
    negfile = corpora / "geneontology-2004.txt"

#posfile = corpora / "genedrug-small.txt"
