"""Configuration for MedScanner

@author: Graham Poulter
                                   

Import this module to set up library paths and configure locations to
files used by the scanner, as well as configuration options.

These are the most useful options for users:

* base.working
* query.{limit, pseudocount, posfile, outputdb}
* validation.{recall, nfolds, pseudocount, posfile, negfile}

"""

import logging
import sys
from new import module
from path import path

__all__ = [ "base", "medline", "genedrug", "update", "query", "validation" ]

#### Basic configuration used by other modules

b = base = module("config_base")
## Path to source directory
b.src = path(sys.argv[0]).dirname().abspath().parent
## Path to top 
b.top = b.src.parent
## Path to bin
b.bin = b.src / "bin"
## Path to lib
b.lib = b.src / "lib"
## Path to templates
b.templates = b.lib / "templates"
## Path to stylesheet
b.stylesheet = b.templates / "style.css"
## Path to base for working data
b.working = b.top / "data"
#b.working = b.top / "testing"
## Path to output files
b.output = b.working / "output"
## Path to cache directory
b.cache = b.working / "cache"
## Path to web output directory
b.weboutput = b.src / "www" / "htdocs" / "output"
## Path to corpora
b.corpora = b.working / "corpora"

#### MEDLINE DATABASE CONFIGURATION

m = medline = module("config_medline")
## Path to BSD DB environment home directory
m.db_home = b.cache / "db_home"
## Path to DB of article objects
m.articledb = b.cache / "articles.db"
## Path to list of article IDs
m.articlelist = b.cache / "articles.txt"
## Path to DB of term features for each article
m.featuredb = b.cache / "features.db"
## Path to to find MeSH term mapping database
m.meshdb = b.cache / "meshdb.txt"
## Path to term counts pickle
m.termcounts = b.cache / "termcounts.pickle"
## Path to list files already processed
m.processed = b.cache / "processed.txt"

#### GENEDRUG CONFIGURATION

gd = genedrug = module("config_genedrug")
## Path to GAPScore results in
gd.gapscore = b.cache / "gapscore.db"
## Path to gene-drug co-occurrence cache
gd.genedrug = b.cache / "genedrug.db"
## Path to pickled drug table (generated by drugtable.py)
gd.drugtable = b.working / "drug-table" / "drugtable.pickle"

#### UPDATE CONFIGURATION (for update.py)

u = update = module("config_update")
## Path to compressed MEDLINE 
u.medline = b.working / "medline"
## Number of files to process between saving
u.save_delay = 10
## Pickle with MeSH terms to exlude (generated by meshexcludes.py)
u.mesh_excludes = b.working / "mesh-excludes" / "meshexcludes.pickle"
## Pickle with MeSH synonyms (generated by meshexcludes.py)
u.mesh_synonyms = b.working / "mesh-excludes" / "meshsynonyms.pickle"

#### QUERY CONFIGURATION (for query.py)

q = query = module("config_query")
## Path to file of positive training PMIDs (processed by fixcorpora.py)
q.posfile = None
## Integer for maximum number of results (may be fewer due to threshold)
q.limit = 10000
## Float for minimum score threshold
q.threshold = 0
## Per-term pseudocount to use
q.pseudocount = 0.01
## Prefix for result report files
q.prefix = b.output / "results"
## Path to database output
q.outputdb = None
## Stylesheet for report
q.stylesheet = b.stylesheet

#### VALIDATOR CONFIGURATION (for validate.py)

v = validation = module("config_validation")
## Positive training set (processed by fixcorpora.py)
v.posfile = None
## Negative training / testing set (processed by fixcorpora.py)
v.negfile = None
## Validation folds to use
v.nfolds = 5
## Pseudocounts to use
v.pseudocount = 0.01
## Whether to do a gene-drug co-occurrence filter
v.genedrug = False
## Whether to use Daniel's 10^-8 pseudocounts
v.daniel = False
## Prefix to use for report file
v.prefix = b.output / "validation"
## Stylesheet for report
v.stylesheet = b.stylesheet

#### Logging configuration

logging.basicConfig(
    level    = logging.DEBUG,
    datefmt  = "%H:%M:%S",
    format   = "%(asctime)-9s %(levelname)-8s %(message)s",
)

#### Configure Python search path

sys.path.insert( 0, base.lib )

#### Configure which data sets to ues

#dataset = "aids-vs-500k"
#dataset = "cur-vs-500k"
#dataset = "cur-vs-med"
#dataset = "cur-vs-go4"
#dataset = "full-vs-500k"
#dataset = "full-vs-go4"
dataset = "old-vs-go4"
#dataset = "old-vs-500k"

q.prefix = b.output / (dataset+"-result")
v.prefix = b.output / (dataset+"-valid")

if dataset == "aids-vs-500k":
    v.posfile = b.corpora / "ncbi-aids-bioethics-Oct06-1k.txt"
    v.negfile = b.corpora / "medline-500k.txt"
if dataset == "cur-vs-med":
    v.posfile = b.corpora / "pharmgkb-Oct06.txt"
    v.negfile = m.articlelist
if dataset == "cur-vs-500k":
    v.posfile = b.corpora / "pharmgkb-Oct06.txt"
    v.negfile = b.corpora / "medline-500k.txt"
if dataset == "cur-vs-go4":
    v.posfile = b.corpora / "pharmgkb-Oct06.txt"
    v.negfile = b.corpora / "geneontology-2004.txt"
if dataset == "full-vs-500k":
    v.posfile = b.corpora / "pharmgkb-full-Jan06.txt"
    v.negfile = b.corpora / "medline-500k.txt"
if dataset == "full-vs-go4":
    v.posfile = b.corpora / "pharmgkb-full-Jan06.txt"
    v.negfile = b.corpora / "geneontology-2004.txt"
if dataset == "old-vs-go4":
    v.posfile = b.corpora / "pharmgkb-2004.txt"
    v.negfile = b.corpora / "geneontology-2004.txt"
if dataset == "old-vs-500k":
    v.posfile = b.corpora / "pharmgkb-2004.txt"
    v.negfile = b.corpora / "medline-500k.txt"

q.posfile = v.posfile
