.. -*- mode: rst -*-

=======================================================
MScanner: a classifier for retrieving Medline citations
=======================================================

:Author: Graham Poulter
                                

Below are instructions for running a local copy of the MScanner_ web
service.  One might want to install it locally for the purpose of
modifying the source code, or to integrate the code into a natural
language processing pipeline.  

.. _MScanner: http://mscanner.stanford.edu


Installation
============

  * Unzip the source code into some directory, which I will refer to
    as ``INSTALL``.  On a Unix-like system, ``INSTALL`` might be
    ``/home/joe/MScanner``.

  * Create any directories required by the `Directory Layout`_.

  * Add ``INSTALL/mscanner`` to the ``PYTHONPATH`` environment
    variable.  You should now be able to start a python shell and type
    ``import mscanner``.

  * Run GNU ``make`` in the directory ``INSTALL/mscanner/fastscores``,
    to compile _FeatureScores.c and _ScoreCalculator.c, for much
    faster queries.

  * Install the Dependencies_ 

  * Set up the Configuration_ file, ``INSTALL/mscanner/configuration.py``
  
  * Place `Input Data`_ in the correct locations.

  * Run the unit tests by executing ``python __init__.py`` in
    ``INSTALL/mscanner/tests``.

  * Run ``INSTALL/mscanner/scripts/update.py`` to parse Medline and create
    the databases in ``INSTALL/data/cache``.

  * Start the `Web Interface`_, or create custom query and validation
    operations by editing query.py and validate.py in
    ``INSTALL/mscanner/scripts``.  The other scripts in the directory
    have documentation in their source code.


Directory Layout
================

The structure of  the ``INSTALL`` directory after extraction is as follows::
                                          
     mscanner/help -- INSTALL.txt
     mscanner/help/api -- Generated API documentation in HTML
     mscanner/core -- The main modules for query and validation
     mscanner/fastscores -- C programs for accelerated classification
     mscanner/htdocs -- Web interface
     mscanner/htdocs/static -- HTML files and other static data
     mscanner/htdocs/static/output -- Where web interface results are placed
     mscanner/medline -- Modules for indexing Medline records
     mscanner/scripts -- Utilities for performing operations
     mscanner/tests -- Unit tests

Additionally create these directories to hold data::

     data/corpora -- Input files
     data/medline -- The downloaded Medline distribution
     data/queue -- Web interface task queue
    

Input Data
==========

MScanner requires the Medline_ Baseline distribution to generate the
database.  Please contact the National Libraries of Medicine for a
license. MScanner looks for the .xml.gz files of the FTP distribution
of Baseline in ``INSTALL/data/medline`` (the ``rc.medline`` parameter)

Training the classifier requires PubMed IDs of relevant training
examples.  See the MScanner publication for more information.  Data
sets should be placed in ``INSTALL/data/corpora`` (set by the
``rc.corpora`` parameter)

.. _Medline: http://www.nlm.nih.gov/bsd/licensee/2007_stats/baseline_doc.html


Dependencies
============

Python_
 The Python programming language (version 2.5)

Numpy_
 A library for Matlab-like vector calculations.
 
path_
 A Python module that allows one to treat directory paths as objects.

pysqlite_
 Latest version of the Python interface to the SQLite database engine.

Cheetah_
 Templating engine for generating HTML outputs.

`Gnuplot-py`_
 For plotting cross validation output graphs. Gnuplot-py in turn
 requires Gnuplot_ to be installed.

`web.py`_
 Python web application framework (for running the web interface).

.. _Python: http://www.python.org
.. _Numpy: http://numpy.scipy.org
.. _path: http://www.jorendorff.com/articles/python/path
.. _Cheetah: http://www.cheetahtemplate.org
.. _pysqlite: http://pysqlite.org/
.. _Gnuplot-py: http://gnuplot-py.sourceforge.net/
.. _Gnuplot: http://www.gnuplot.info/ 
.. _web.py: http://webpy.org


Web Interface
=============

The web interface requires `web.py`_ to be installed. web.py has a
built-in webserver, but the pages can also behave as CGI scripts.

Start the queue process:

  Run ``INSTALL/mscanner/htdocs/queue.py`` to start a long-running
  process that watches ``INSTALL/data/queue`` for descriptor files
  deposited by the web interface.  When a descriptor file arrives, the
  queue.py performs the requested query/validation task, and deposits
  the output in the static HTML directory.

Start the internal web server:

  Execute ``INSTALL/mscanner/htdocs/controller.py``.  The web server
  will run on ``http://localhost:8080``.  The web interface takes user
  input, reports on status of MScanner, and places descriptor files in
  the queue directory on submission of a query.

Or: run the web interface under Apache:

   Configure Apache to use ``controller.py`` as a CGI script.  Without
   URL rewriting you will have to access pages as 
   ``http://localhost/controller.py/form``.


Re-generating API Documentation
-------------------------------

API documentation is linked from the MScanner_ home page, and may
already be provided under ``INSTALL/mscanner/doc/api``.

To re-generate the documentation, first install Epydoc_ , then navigate to
the the ``INSTALL/mscanner`` directory, and run::

    epydoc.py --config=doc/epydoc.txt

To convert ``INSTALL/mscanner/doc/INSTALL.txt`` (which is written in the
reStructuredText_ format) into ``INSTALL.html``, first install
Docutils_, then navigate to ``INSTALL/mscanner/doc`` and run::

    rst2html.py INSTALL.txt > INSTALL.html

.. _Epydoc: http://epydoc.sourceforge.net
.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Docutils: http://docutils.sourceforge.net
  

Configuration
=============

Configuration is in ``INSTALL/mscanner/configuration.py``, which has
MScanner resource parameters (and their descriptions). The ``rc``
object holds all the global constants used by MScanner.  

The most important options to set rc.articles_path (path
for storing the article database) and rc.medline_path (where to find
the Medline XML files).

See the scripts in ``INSTALL/mscanner/scripts`` for examples of using
the MScanner API.  Help on the API is in
``INSTALL/mscanner/help/api``.

