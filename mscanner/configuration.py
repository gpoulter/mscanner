"""Configuration for MedScanner

Configuration options are stored in the "rc" variable in this module, which is a
utils.RCStorage object.

@note: Some RC parameters, like directory paths, provide configuration
for things which tend not to change between invocations and data sets.
Others, like pseudocount or path to the input file, are an easy
way to provide centralised defaults for those options.  Instantiated
objects, like FeatureInfo, keep their own pseudocount etc.
in the instance (passed via __init__), which may differ from the RC one.

@note: For parameters which depend on others, prefix with lambda:
so that the dependency is updated dynamically.
"""

                                               
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

from path import path
from mscanner.support.storage import RCStorage

rc = RCStorage()

#### PATH CONFIGURATION

## Path to MScanner source directory
rc.sources = path(__file__).dirname().parent
## Path to directory for report templates
rc.templates = lambda: rc.sources / "mscanner" / "templates"
## Path to working information directory
rc.working = lambda: rc.sources.parent / "data"
## Path to cache directory
rc.cache = lambda: rc.working / "cache"
## Path to web directory
rc.htdocs = lambda: rc.sources / "htdocs"
## DB environment directory
rc.db_env_home = lambda: rc.cache / "db_home"
## DB of article objects
rc.articledb = lambda: rc.cache / "articles.db"
## Path for list of of article IDs
rc.articlelist = lambda: rc.cache / "articles.txt"
## Path for number of articles
rc.narticles = lambda: rc.cache / "narticles.txt"
## Path for DB of term features for each article
rc.featuredb = lambda: rc.cache / "features.db"
## Path for binary stream of PMIDs and feature arrays
rc.featurestream = lambda: rc.cache / "features.stream"
## Path for feature<->ID mapping
rc.featuremap = lambda: rc.cache / "featuremap.txt"
## Path for list of files already processed
rc.processed = lambda: rc.cache / "processed.txt"
## Path for directory with compressed Medline
rc.medline = lambda: rc.working / "medline"
## Path for log file
rc.logfile = lambda: rc.cache / "lastlog.txt"
## Path for e-mail alert file
rc.emails_path = lambda: rc.cache / "emails.txt"
## Path for stylesheet for reports
rc.stylesheet = lambda: rc.templates / "style.css"
## Path to corpora directory
rc.corpora = lambda: rc.working / "corpora"
## Path to outputs directory for queries
rc.query_report_dir = lambda: rc.working / "query" / rc.dataset
## Path to outputs directory for validation
rc.valid_report_dir = lambda: rc.working / "validation" / rc.dataset
## Path to outputs directory for web submissions
rc.web_report_dir = lambda: rc.htdocs / "static" / "output"
## Path to the cscore program
rc.cscore_path = lambda: rc.sources / "mscanner" / "cscore" / "cscore"
## Path to the descriptor queue
rc.queue_path = lambda: rc.working / "queue"


### GENE-DRUG DETECTION
rc.genedrug = lambda: rc.working / "genedrug"
## Path to GAPScore results db
rc.gapscore_db = lambda: rc.genedrug / "gapscore.db"
## Path to gene-drug co-occurrence cache
rc.genedrug_db = lambda: rc.genedrug / "genedrug.db"
## Path to drug table
rc.drugtable = lambda: rc.genedrug / "drugtable.txt"
## Path to CSV table of co-occurrences
rc.genedrug_csv = lambda: rc.genedrug / "pharmdemo.csv"
## Path to SQL table of co-occurrences
rc.genedrug_sql = lambda: rc.genedrug / "pharmdemo.sql"

### ALL OUTPUTS
## Name of index file
rc.report_index = path("index.html")
## Name of descriptor file
rc.report_descriptor = path("descriptor.txt")
### QUERY
## Name of list of scores of input PMIDs
rc.report_input_scores = path("input_scores.txt")
## Name of list of broken input PMIDs
rc.report_input_broken = path("input_broken.txt")
## Name of result score file
rc.report_result_scores = path("result_scores.txt")
## Name of file with citation records for the input
rc.report_input_citations = path("input_citations.html")
## Name of file with citation records for the output
rc.report_result_citations = path("result_citations.html")
### RETRIEVAL TEST
## Name of file with list of testing PMIDs
rc.report_retrieval_test_pmids = path("retrieval_test.txt")
## Name of file with cumulative count of retrieved test PMIDs
rc.report_retrieval_stats = path("retrieval_stats.txt")
## Name of retrieval vs rank graph
rc.report_retrieval_graph = path("retrieval.png")
### VALIDATION
## Name of file with positive PMIDs and scores
rc.report_positives = path("positives.txt")
## Name of file with broken positive PMIDs
rc.report_positives_broken = path("positives_broken.txt")
## Name of file with negative PMIDs and scores
rc.report_negatives = path("negatives.txt")
## Name of file with excluded negative PMIDs
rc.report_negatives_exclude = path("negatives_exclude.txt")
## Name of term score file
rc.report_term_scores = path("term_scores.csv")
## Name of histogram file
rc.report_artscores_img = path("artscores.png")
## Name of feature score density file
rc.report_featscores_img = path("featscores.png")
## Name of ROC curve file
rc.report_roc_img = path("roc.png")
## Name of PR curve file
rc.report_prcurve_img = path("prcurve.png")
## Name of PRF vs threshold file
rc.report_fmeasure_img = path("fmeasure.png")

#### Non-Path parameters

## Whether to use transactions while updating
rc.use_transactions = False
## Number of seconds to pause before next file while updating
rc.save_delay = 6
## Server for sending e-mails
rc.smtpserver = "smtp.uct.ac.za" # "smtp.stanford.edu"
## Proportion of data to use in retrieval test
rc.retrieval_test_prop = 0.1
## Timestamp to use in the output
rc.timestamp = None
## Email to send website queries to
rc.webmaster_email = "xxxxxxxxxxxxxxxxxxxxxxxx"
## Base directory for the website
rc.siteurl = "http://mscanner.stanford.edu"
## If true, link to .js and .css files instead of including them
rc.linkHeaders = False;

## Parameters affecting FeatureInfo 

## Per-term pseudocount to use (None for background frequency)
rc.pseudocount = None
## Types of features to exclude
rc.exclude_types = []
## Method name for calculating feature probabilities
rc.getFrequencies = "getProbabilitiesBayes"
## Method name for calculating mask after scores (may be None)
rc.getPostMask = None

#### Web-configurable parameters

## Name of dataset
rc.dataset = "default"
## Float for minimum score threshold
rc.threshold = 10
## Integer for maximum number of results (may be fewer due to threshold)
rc.limit = 1000
## Number of validation folds to use
rc.nfolds = 10
## 0<Alpha<1.  Alpha=0.5 maximises standard F-Measure.
rc.alpha = 0.5
## Number of negatives to use
rc.numnegs = 500000
## Number of citations per output file
rc.citations_per_file = 200

## Logging configuration
def initLogger(console=True, logfile=None):
    """Configure logging for MScanner
    
    @param console: If True, log to the console.
    
    @param logfile: File for logging, defaulting to rc.logfile. If False, do
    not log to file. """
    import logging
    # Configure root logger
    rootlog = logging.getLogger()
    rootlog.setLevel(0)
    format = logging.Formatter("%(asctime)-9s %(levelname)-8s %(message)s", "%H:%M:%S")
    # Configure file logging
    if logfile != False:
        filelog = logging.FileHandler(logfile if logfile else rc.logfile, "w")
        filelog.setFormatter(format)
        rootlog.addHandler(filelog)
    # Configure console logging
    if console:
        console = logging.StreamHandler()
        console.setFormatter(format)
        rootlog.addHandler(console)
