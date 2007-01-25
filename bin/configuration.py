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
sys.path.insert(0, lib)
## Path to templates
templates = lib / "templates"
## Path to base for working data
working = top / "data"
#working = top / "testing"
## Path to cache directory
cache = working / "cache"

#### DATABASE CONFIGURATION

## Path to stylesheet for query/validation reports
stylesheet = templates / "style.css"
## Path to BSD DB environment home directory
db_home = cache / "db_home"
## Path to DB of article objects
articledb = cache / "articles.db"
## Path to list of article IDs
articlelist = cache / "articles.txt"
## Path to DB of term features for each article
featuredb = cache / "features.db"
## Path to feature<->ID mapping
featuremap = cache / "featuremap.txt"
## Path to list files already processed
processed = cache / "processed.txt"
## Whether to use transactions while updating
use_transactions = False
## Path to compressed MEDLINE 
medline = working / "medline"
## Number of files to process between saving
save_delay = 5

#### GENEDRUG CONFIGURATION

## Path to GAPScore results in
gapscore = working / "genedrug" / "gapscore.db"
## Path to gene-drug co-occurrence cache
genedrug = working / "genedrug" / "genedrug.db"
## Path to pickled drug table (generated by drugtable.py)
drugtable = working / "genedrug" / "drugtable.txt"

#### Shared configuration (by query & validation) 

## Per-term pseudocount to use
pseudocount = 0.01
## Path to output files
output = working / "output"
#output = src / "www" / "htdocs" / "output"
## Path to corpora
corpora = working / "corpora"
## Prefix for report files (output/dataset)
reportdir = None
## File name for positive documents
posfile = "positives.txt"
## Types of features to exclude
exclude_types = None # ["issn"]
## Path to file with progress statistics
statfile = cache / "mscanner.pid"
## File of e-mail addresses to alert on completion
mailer = cache / "emails.txt"
## Server for sending e-mails
smtp_server = "smtp.uct.ac.za" # "smtp.stanford.edu"
## File name for report index page
index_file = "index.html"
## Base file name for term scores
term_scores_name = "term_scores.csv"

#### QUERY CONFIGURATION (for query.py)

## Integer for maximum number of results (may be fewer due to threshold)
limit = 10000
## Float for minimum score threshold
threshold = 0
## Path to database output
outputdb = None
## Base file name for results
query_results_name = "result_scores.txt"
## Base file name for input citations
input_citations = "input_citations.html"
## Base file name for result citations
result_citations= "result_citations.html"

#### VALIDATOR CONFIGURATION (for validate.py)

## Validation folds to use
nfolds = 0
## 0<Alpha<1.  Alpha=0.5 maximises standard F-Measure.
alpha = 0.5
## File name for negative documents
negfile = "negatives.txt"
## Maximum number of negatives
numnegs = 500000
## Whether to do a gene-drug co-occurrence filter
dogenedrug = False
## Whether to use Daniel's 10^-8 pseudocounts
dodaniel = False
## File name for histogram
hist_img = "histogram.png"
## File name for ROC curve
roc_img = "roc.png"
## File name for PR curve
p_vs_r_img = "prcurve.png"
## File name for PRF vs threshold
pr_vs_score_img = "prscore.png"

#### DATA SET CONFIGURATION

def configure_query(dataset, pseudocount, limit, threshold, pos):
    import configuration as c
    c.dataset = dataset
    c.pseudocount = pseudocount
    c.limit = limit
    c.threshold = thresholdp
    c.reportdir = output / dataset
    pos.copy(reportdir/posfile)

def configure_validation(dataset, pseudocount, numnegs, alpha, pos, neg):
    c.dataset = datset
    c.pseudocount = pseudocount
    c.numnegs = numnegs
    c.alpha = alpha
    c.reportdir = output / dataset
    pos.copy(posfile)
    if isinstance(neg, basestring):
        neg.copy(reportdir/negfile)
    elif isinstance(neg, int):
        (reportdir/negfile).write_lines(random.sample(articlelist.lines(), numnegs))
    
def choose_query(dataset):
    import configuration as c
    c.dataset = dataset
    c.reportdir = output / (dataset+"-query")
    if dataset == "pg04":
        pos = "pharmgkb-2004.txt"
    if dataset == "pg06":
        pos = "pharmgkb-Oct06.txt"
    if dataset == "aids":
        pos = "aids-bioethics-Oct06.txt"
    if dataset == "radiology":
        pos = "daniel-radiology.txt"
    if dataset == "mscanner":
        pos = "mscanner-bibliography.txt"
    if dataset == "gdsmall":
        pos = "genedrug-small.txt"
    if not isinstance(pos, path):
        pos = corpora / pos
    if not reportdir.isdir():
        reportdir.mkdir()
        pos.copy(reportdir/posfile)
    
def choose_validation(dataset):
    import configuration as c
    c.dataset = dataset
    c.reportdir = output / (dataset+"-valid")
    # Comparing to daniel's old method
    if dataset == "pg04-vs-30k":
        pos = "pharmgkb-2004.txt"
        neg = "medline06-30k.txt"
    if dataset == "pg04-vs-30k-dan":
        pos = "pharmgkb-2004.txt"
        neg = "medline06-30k.txt"
        c.exclude_feats = ["issn"]
        c.dodaniel = True
    if dataset == "pg06-vs-500k-dan":
        pos = "pharmgkb-Oct06.txt"
        neg = "medline06-500k.txt"
        c.exclude_feats = ["issn"]
        c.dodaniel = True
    if dataset == "pg06-vs-500k-noissn":
        pos = "pharmgkb-Oct06.txt"
        neg = "medline06-500k.txt"
        c.exclude_feats = ["issn"]
    # Primary results
    if dataset == "aids-vs-500k":
        pos = "aids-bioethics-Oct06.txt"
        neg = "medline06-500k.txt"
    if dataset == "pg06-vs-500k":
        pos = "pharmgkb-Oct06.txt"
        neg = "medline06-500k.txt"
    if dataset == "radiology-vs-500k":
        pos = "daniel-radiology.txt"
        neg = "medline06-500k.txt"
    if dataset == "random10k-vs-500k":
        pos = "random10k-06.txt"
        neg = "medline06-500k.txt"
    # Other experiments
    if dataset == "mscanner-vs-500k":
        pos = "mscanner-bibliography.txt"
        neg = "medline06-500k.txt"
    if dataset == "pg06-vs-med06":
        pos = "pharmgkb-Oct06.txt"
        neg = articlelist
    if dataset == "gdsmall-vs-med06":
        pos = "genedrug-small.txt"
        neg = articlelist
    if dataset == "aids1k-vs-500k":
        pos = "aids-bioethics-Oct06-1k.txt"
        neg = "medline06-500k.txt"
    if dataset == "pg04-vs-go4":
        pos = "pharmgkb-2004.txt"
        neg = "geneontology-2004.txt"
    if dataset == "pg04-vs-500k":
        pos = "pharmgkb-2004.txt"
        neg = "medline06-500k.txt"
    if dataset == "pgdan-vs-500k":
        pos = "pharmgkb-daniel.txt"
        neg = "medline06-500k.txt"
    if not isinstance(pos, path):
        pos = corpora / pos
    if not isinstance(neg, path):
        neg = corpora / neg
    if not reportdir.isdir():
        reportdir.mkdir()
        pos.copy(reportdir/posfile)
        neg.copy(reportdir/negfile)
