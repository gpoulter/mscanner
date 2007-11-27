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

#### BASIC MSCANNER PATHS

## Path to MScanner source directory
rc.sources = path(__file__).dirname()
## Path to web directory
rc.htdocs = lambda: rc.sources / "htdocs"
## Path to directory for report templates
rc.templates = lambda: rc.sources / "core" / "templates"
## Path to list of stop words
rc.stopwords = lambda: rc.sources / "stopwords.txt"
## Path to working information directory
rc.working = lambda: rc.sources.parent / "data"

### WORKING DIRECTORY PATHS

## Path for directory with compressed Medline
rc.medline = lambda: rc.working / "medline"
## Path to corpora directory
rc.corpora = lambda: rc.working / "corpora"
## Path to outputs directory
rc.web_report_dir = lambda: rc.htdocs / "static" / "output"
## Path to the descriptor queue
rc.queue_path = lambda: rc.working / "queue"

## CACHE DIRECTORY FILES

## Path to cache directory
rc.cache = lambda: rc.working / "cache"
## Path to database environment directory
rc.db_env_home = lambda: rc.cache / "db_home"
## Path to Article object database
rc.articledb = lambda: rc.cache / "articles.db"
## Path to list of of article IDs
rc.articlelist = lambda: rc.cache / "articles.txt"
## Path to track Medline files that have already been added
rc.tracker = lambda: rc.cache / "processed.txt"
## Path for log file
rc.logfile = lambda: rc.cache / "lastlog.txt"

## FEATURE DATABASE FILES (16 and 32-bit)

## Path for DB of MeSH-features for each article
rc.featuredb_mesh = lambda: rc.cache / "features_mesh.db"
## Path for DB of all-features for each article
rc.featuredb_all = lambda: rc.cache / "features_all.db"
## Path for binary stream of PMIDs and MeSH-feature arrays
rc.featurestream_mesh = lambda: rc.cache / "features_mesh.stream"
## Path for binary stream of PMIDs and All-feature arrays
rc.featurestream_all =  lambda: rc.cache / "features_all.stream"
## Path for feature<->ID mapping for MeSH features
rc.featuremap_mesh = lambda: rc.cache / "featuremap_mesh.txt"
## Path for feature<->ID mapping for All features
rc.featuremap_all =  lambda: rc.cache / "featuremap_all.txt"


### ALL REPORTS

## Name of index file
rc.report_index = path("index.html")
## Name of descriptor file
rc.report_descriptor = path("descriptor.txt")
## Name of term score file
rc.report_term_scores = path("terms.csv")
## Name of file to echo the log to
rc.report_logfile = path("logging.txt")

### QUERY OUTPUTS

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
rc.report_result_all_zip = rc.report_result_all + ".zip"

### VALIDATION OUTPUTS

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

## Parameters affecting FeatureScores 

## Per-term pseudocount to use (None for background frequency)
rc.pseudocount = None
## Types of features to exclude
rc.exclude_types = []
## Method name for calculating feature probabilities
rc.make_scores = "scores_bayes"
## Method name for calculating mask after scores (may be None)
rc.get_postmask = None