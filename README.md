# MScanner

MScanner retrieves Medline abstracts using machine learning. The input is a list of
PubMed IDs of articles on some topic, and the output is a page of articles most relevant 
to the topic as inferred from the examples. I wrote MScanner for my Masters 
Thesis on Rapid statistical classification on the Medline database of biomedical literature.

It uses some search engine techniques to speed up classification, ranking about 500k
abstracts per second single-threaded.

* The service is hosted at http://mscanner.stanford.edu
* The issue tracker is on [Trello](https://trello.com/board/mscanner-board/4f0d30714ac4a64f673a2dda)

Automatically exported from https://code.google.com/p/mscanner/
