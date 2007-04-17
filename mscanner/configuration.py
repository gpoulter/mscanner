"""Configuration for MedScanner

Import this module for configuration options.

                                   

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

import logging
from path import path
import time

#### Logging configuration

logging.basicConfig(
    level    = logging.DEBUG,
    datefmt  = "%H:%M:%S",
    format   = "%(asctime)-9s %(levelname)-8s %(message)s",
)

#### TOP-LEVEL CONFIGURATION 
## (other options can reference these but no others)

## Path to MScanner source directory
src = path(__file__).dirname().parent
## Path to top MScanner directory
top = src.parent
## Path to directory for report templates
templates = src / "mscanner" / "templates"
## Path to working information directory
working = top / "data" # "testing"
## Path to cache directory
cache = working / "cache"
## Path to web directory
htdocs = src / "htdocs"

#### DATABASE CONFIGURATION

## Path to stylesheet for reports
stylesheet = templates / "style.css"
## Path to DB environment directory
db_home = cache / "db_home"
## Path to DB of article objects
articledb = cache / "articles.db"
## Path to list of article IDs
articlelist = cache / "articles.txt"
## Path to count of number of articles
narticles = cache / "narticles.txt"
## Path to DB of term features for each article
featuredb = cache / "features.db"
## Path to binary stream of PMIDs and feature arrays
featurestream = cache / "features.stream"
## Path to feature<->ID mapping
featuremap = cache / "featuremap.txt"
## Path to list of files already processed
processed = cache / "processed.txt"
## Whether to use transactions while updating
use_transactions = False
## Path to directory with compressed Medline
medline = working / "medline"
## Number of files to process between saving
save_delay = 2

#### GENEDRUG CONFIGURATION

## Path to GAPScore results db
gapscore = working / "genedrug" / "gapscore.db"
## Path to gene-drug co-occurrence cache
genedrug = working / "genedrug" / "genedrug.db"
## Path to pickled drug table (generated by drugtable.py)
drugtable = working / "genedrug" / "drugtable.txt"

#### Shared configuration (by query & validation) 

## Per-term pseudocount to use
pseudocount = 0.01
## Whether to set to zero negative-only features with positive scores 
cutoff = False
## Path to corpora directory
corpora = working / "corpora"
## Prefix for report files (output/dataset)
reportdir = None
## Types of features to exclude
exclude_types = []
## Time at which this file was imported
timestamp = time.time()
## Path to status file
statfile = cache / "mscanner.pid"
## Path to e-mail alert file
emails_path = cache / "emails.txt"
## Server for sending e-mails
smtp_server = "smtp.uct.ac.za" # "smtp.stanford.edu"
## Name of file with positive PMIDs
posname = "positives.txt"
## Name of index file
index_html = "index.html"
## Name of term score file
term_scores = "term_scores.csv"

#### QUERY CONFIGURATION (for query.py)

## Integer for maximum number of results (may be fewer due to threshold)
limit = 10000
## Float for minimum score threshold
threshold = 0
## Path to database output
outputdb = None
## Name of result score file
result_scores = "result_scores.txt"
## Name of file with input citations
input_citations = "input_citations.html"
## Name of file with result citations
result_citations = "result_citations.html"
## Path to outputs directory
query_output = working / "query"

#### VALIDATOR CONFIGURATION (for validate.py)

## Validation folds to use
nfolds = 10
## 0<Alpha<1.  Alpha=0.5 maximises standard F-Measure.
alpha = 0.5
## Maximum number of negatives
numnegs = 500000
## Whether to do a gene-drug co-occurrence filter
dogenedrug = False
## Whether to use Daniel's 10^-8 pseudocounts
dodaniel = False
## Name of file with negative PMIDs
negname = "negatives.txt"
## Name of histogram file
artscores_img = "artscores.png"
## Name of feature score density file
featscores_img = "featscores.png"
## Name of ROC curve file
roc_img = "roc.png"
## Name of PR curve file
prcurve_img = "prcurve.png"
## Name of PRF vs threshold file
fmeasure_img = "fmeasure.png"
## Path to outputs directory
valid_output = working / "validation"
