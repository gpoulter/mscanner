"""Configuration for MedScanner

Some parameters change while the program is running, such as the query and
validation parameters.

Parameters can be made to depend on others by using lambda, which the rc object
auto-calls those so they appear as normal data attributes."""


from path import path
from mscanner.core.Storage import RCStorage


                                     
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


rc = RCStorage()
"""Global configuration options"""

## DYNAMIC PATHS

## This file is at the top source directory
rc.sources = path(__file__).dirname()
## Path to root directory containing sources and data
rc.root = lambda: rc.sources.parent
## Path to data directory (inputs and databases)
rc.working = lambda: rc.root / "data"
## Path to ArticleData (with FeatureData subdirs)
rc.articles_home = lambda: rc.working / "articles_bmc"
## Path to report templates
rc.templates = lambda: rc.sources / "core" / "templates"
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

## ARTICLE DATABASE FILES (in articles_home)

## Path to database environment directory
rc.db_env_home = path("db_home")
## Path to Article object database
rc.articledb = path("articles.db")
## Path to list of of article IDs
rc.articlelist = path("articles.txt")
## Path to track Medline files that have already been added
rc.tracker = path("processed.txt")
## Path for log file
rc.logfile = path("lastlog.txt")

## FEATURE DATABASE FILES (in a subdir of articles_home)

## Base name for DB of MeSH-features for each article
rc.featuredb = path("features.db")
## Base name for binary stream of PMIDs and MeSH-feature arrays
rc.featurestream = path("features.stream")
## Base name for feature<->ID mapping for MeSH features
rc.featuremap = path("featuremap.txt")

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

## Name of list of scores of input PMIDs
rc.report_input_scores = path("inputs.txt")
## Name of list of broken input PMIDs
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

## Name of file with positive PMIDs and scores
rc.report_positives = path("positives.txt")
## Name of file with broken positive PMIDs
rc.report_positives_broken = path("positives_broken.txt")
## Name of file with negative PMIDs and scores
rc.report_negatives = path("negatives.txt")
## Name of file with excluded negative PMIDs
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
## Name of predicted precision graph
rc.report_prediction_img = path("prediction.png")

#### Non-Path parameters

## Server for sending e-mails
rc.smtpserver = "smtp.uct.ac.za" # "smtp.stanford.edu"
## Email to send website queries to
rc.webmaster_email = "xxxxxxxxxxxxxxxxxxxxxxxx"
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
rc.randseed = 124

## Parameters affecting FeatureScores 

## Exclude features occurring fewer times in the data
rc.mincount = 1
## Exclude features of theses classes. e.g. ["mesh","issn"]
rc.class_mask = []
## Exclude features with less than this Information Gain (IG)
rc.min_infogain = 0
## If True, exclude features not occurring in any relevant articles
rc.positives_only = False
## Method name for calculating feature probabilities
rc.scoremethod = "scores_laplace_split"
## Print at most this many feature scores
rc.max_output_features = 50000

## Initialise article stopwords

from mscanner.medline import Article
Article.init_stopwords(rc.stopwords)