I wrote MScanner while working towards my Master of Science degree in bioinformatics.  My [website page](https://sites.google.com/site/grahampoulter/m-sc-thesis) and the  [thesis itself](https://sites.google.com/site/grahampoulter/m-sc-thesis/masters-thesis.pdf?attredirects=0) (PDF) have more information about the research

# Functionality #

MScanner carries out the following functions:
  * Indexes Medline (16 millions documents) into a compact database of features.
  * Trains a Naive Bayes classifier from positive examples, using rest of Medline counts as an approximation to the negative class.
  * Classifies the compressed documents
  * Serves a web application for submitting jobs and viewing results

# Algorithms #

The machine learning is sped up using techniques from  information retrieval, namely a compressed index of numerical features, and a transformation of the Naive Bayes formula for fast evaluation of the classifier score.  In production, MScanner classifies 100,000 documents per second on a single core.

MScanner could be further sped up by caching features persistently in memory, and distributing the classification over multiple cores in a map-reduce fashion.