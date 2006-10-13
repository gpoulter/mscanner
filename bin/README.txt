==========
MedScanner
==========

:Author: Graham Poulter
                                     

Files in bin/
-------------

configuration.py -- Provide options and paths for update/query/validate
drugtable.py -- Parse drug table into a Pickle
meshexludes.py -- Parse MeSH descriptors for exclusions and synonyms
fixcorpora.py -- Make sure positives, negatives are unique, disjoint, non-bad
geneontology.py -- Extract PMIDs from GeneOntology annotations
query.py -- Query MedScanner for results 
update.py -- Add XML directory (or a Pickle) of Articles to database
validate.py -- Do cross-validation and output performance stats
