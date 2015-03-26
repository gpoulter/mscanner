MScanner currently depends on 4-years old versions of these libraries on Python 2.5.  I am in the process of updating MScanner for Python 2.6 and the latest stable versions.

# Dependencies for Ubuntu #

Use `sudo apt-get install <package>` to install the dependencies on Ubuntu:

  * git
  * gnuplot
  * python2.6
  * python2.6-dev
  * python-bsddb3
  * python-cheetah
  * python-gnuplot
  * python-lxml
  * python-numpy
  * python-matplotlib
  * python-webpy
  * sqlite3

# Dependencies for Solaris #

To install MScanner on Solaris will require the following dependencies, listed here from the [OpenCSW](http://www.opencsw.org/) package repository, which has the most up-to-date open-source packages.

  * berkelydb48 (CSWbdb48) -- key-value database for articles
  * berkelydb48\_devel  (CSWbdb48devel) -- to compile bsddb3 4.8.3 for python
  * git (CSWgit) -- source control
  * gnuplot (CSWgnuplot) -- drawing graphs
  * python (CSWpython) -- python 2.6.2
  * python\_dev (CSWpython-dev) -- compiling modules
  * py\_cheetah (CSWpy-cheetah) -- templating engine
  * py\_lxml (CSWpy-lxml) -- xml parsing
  * py\_numpy (CSWpy-numpy) -- array handling
  * py\_webpy (CSWpy-web) -- web framework

The following packages then need to be installed or compiled from source:

  * bsddb3 python package for Berkely DB
  * gnuplot-py python package for interacting with Gnuplot

Note: we use gnuplot-py on Solaris Sparc because Matplotlib is not packaged and I could not get it to compile despite much effort.

The following packages are useful around MScanner or as compilation dependencies:

  * sqlite3 (CSWsqlite3)
  * libreadline6 (CSWlibreadline6)
  * libreadline6\_dev (CSWlibreadline-dev)
  * ack (CSWack) -- handy super-grep
  * tree (CSWtree) -- tree directory listing
  * vim (CSWvim) -- newer vim