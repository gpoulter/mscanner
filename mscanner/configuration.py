"""Configuration for MedScanner

Some parameters change while the program is running, such as the query and
validation parameters.

Parameters can be made to depend on others by using lambda, which the rc object
auto-calls those so they appear as normal data attributes."""


from path import path
from mscanner.core.Storage import RCStorage


rc = RCStorage()
"""Global configuration options"""

## DYNAMIC PATHS

## This file is at the top source directory
rc.sources = path(__file__).dirname().dirname()
## Path to root directory containing sources and data
rc.root = lambda: rc.sources.parent
## Path to data directory (inputs and databases)
rc.working = lambda: rc.root / "data"
## Path to ArticleData (with FeatureData subdirs)
rc.articles_home = lambda: rc.working / "articles_med12"
## Path to report templates
rc.templates = lambda: rc.sources / "mscanner" / "core" / "templates"
## Path to list of stop words
rc.stopwords = lambda: rc.sources / "stopwords.txt"
## Path to Medline XML
rc.medline = lambda: rc.working / "medline"
## Path to corpora
rc.corpora = lambda: rc.working / "corpora"
## Path to web output directory
rc.web_report_dir = lambda: rc.sources / "htdocs" / "static" / "output"
## Path to web descriptor queue
rc.queue_path = lambda: rc.working / "queue"
## URL for the web root followed by a "/"
rc.web_root = "http://mscanner.stanford.edu/"
## File with process ID of the queue
rc.queue_pid = rc.working / "queue_pid"

## ARTICLE DATABASE FILES (in articles_home)

## Path to Article object database
rc.articledb = path("articles.db")
## Path to track Medline files that have already been added
rc.tracker = path("processed.txt")
## Path for log file
rc.logfile = path("lastlog.txt")

## FEATURE DATABASE FILES (in a subdir of articles_home)

## Base name for feature<->ID mapping for MeSH features
rc.featuremap = path("featuremap.sqlite")
## Base name for DB of MeSH-features for each article
rc.featuredb = path("featvectors.sqlite")
## Base name for binary stream of PubMed IDs and MeSH-feature arrays
rc.featurestream = path("features.stream")

### COMMON REPORT FILES

## Name of index file
rc.report_index = path("index.html")
## Name of descriptor file
rc.report_descriptor = path("descriptor.txt")
## Name of term score file
rc.report_term_scores = path("terms.csv")
## Name of file to echo the log to
rc.report_logfile = path("logging.txt")

### QUERY REPORT FILES

## Name of list of scores of input PubMed IDs
rc.report_input_scores = path("inputs.txt")
## Name of list of broken input PubMed IDs
rc.report_input_broken = path("broken.txt")
## Name of result score file
rc.report_result_scores = path("results.txt")
## Name of file with citation records for the input
rc.report_input_citations = path("inputs.html")
## Name of first page with citation records for the output
rc.report_result_citations = path("results.html")
## Name of file with *ALL* outupt citation records
rc.report_result_all = path("all_results.html")
## Name of zip file with ALL output citations records
rc.report_result_all_zip = path(rc.report_result_all + ".zip")

### VALIDATION REPORT FILES

## Name of file with positive PubMed IDs and scores
rc.report_positives = path("positives.txt")
## Name of file with broken positive PubMed IDs
rc.report_positives_broken = path("positives_broken.txt")
## Name of file with negative PubMed IDs and scores
rc.report_negatives = path("negatives.txt")
## Name of file with excluded negative PubMed IDs
rc.report_negatives_exclude = path("negatives_exclude.txt")
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

## Server for sending e-mails
rc.smtpserver = "smtp.stanford.edu"
## Email to send website queries to
rc.webmaster_email = "graham.poulter@gmail.com"
## Base directory for the website
rc.siteurl = "http://mscanner.stanford.edu"
## Whether to link to .js and .css files instead of including them
rc.link_headers = False
## 0<Alpha<1.  Alpha=0.5 maximises standard F-Measure.
rc.alpha = 0.5
## The utility of retrieving a positive article (defaults to N/P)
rc.utility_r = None
## Number of citations per output file
rc.citations_per_file = 250
## Random seed to use for cross validation shuffle (to get the same
## shuffle each time).  Set to None to get a different seed on each run.
#rc.randseed = 124
rc.randseed = None

## Parameters affecting FeatureScores

## Select features with at least this many occurrences
rc.mincount = 0
## Don't select features of theses classes. e.g. ["mesh","issn"]
rc.type_mask = []
## Select features with at least this much Information Gain (IG)
rc.min_infogain = 0
## If True, select only features found in relevant articles
rc.positives_only = False
## Method name for calculating feature probabilities
rc.scoremethod = "scores_laplace_split"
## Print at most this many feature scores
rc.max_output_features = 50000

## Initialise article stopwords

from mscanner.medline import Article
Article.init_stopwords(rc.stopwords)
