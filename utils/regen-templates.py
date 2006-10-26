#!/usr/bin/env python

"""Regenerate templates used for result and validation output

Usage: python regenerate-templates.py

"""

import sys
from subprocess import call
from path import path

top = path(sys.argv[0]).dirname().abspath().parent
python = path("/usr/bin/python")
templates = top/"lib"/"templates"
preppy = templates/"preppy.py"

call([python,preppy,"compile",templates])
