# Installing #

Clone the git repository (use the "develop" branch) from the Source tab, and [install the dependencies](InstallingDependencies.md), ensuring that you clone it to a directory named "mscanner"

Download the [NLM Sample Data](http://www.nlm.nih.gov/bsd/sample_records_avail.html), which is useful for testing.

For full functionality, apply for a free  [Medline license](http://www.nlm.nih.gov/databases/license/license.html) from the NLM and download the 12GB zipped Medline download from the NLM FTP servers.  Note that the license does not grant any rights to re-distribute downloaded data in bulk format.

Ensure gcc and libc6-dev are available, and run `make` in the `fastcores` directory to build the native `_FeatureCounter` and `_ScoreCalculator` executables.

Ensure "cheetah" executable is in PATH, and run `make` in the `htdocs/templates` directory to compile the Cheetah .tmpl templates to Python files.

If epydoc is installed, run `epydoc --config=help/epydoc.txt` to auto-generate API documentation.

# Importing Data #


From the root of the checkout, do the following:

  * Do `export PYTHONPATH=$(abspath ..)`
  * Place downloaded `*.gz` Medline files in `../data/medline/`
  * Create `../data/articles_med12/`
  * Index the files by running `python scripts/update.py medline`

# Caveats #

MScanner was written with libraries circa 2008 and for Python 2.5, so its a wonder it works as far as importing documents.

I've made some tweaks to work with Python 2.6/2.7 as well (still 2.5-compatible) and to work with newer SQLite.  It uses web.py 0.2 for the web interface and does not work with web.py 0.3+ (I'll fix this once I upgrade web.py on live).  Other dependencies seem to work fine.

The configuration system is convoluted, so its better to place files as described than to edit the settings in configuration.py.

It lacks setup.py and lacks an "mscanner" sub-folder, so the checkout directory itself must be named "mscanner" and the parent of the checkout must be manually added to PYTHONPATH.