"""Manage a status file, which posts program status to any listeners

Format is as follows:

0: PID
1: Progress
2: Total
3: Dataset Identifie
4: Timestamp (via time.time())

@var owner: True if we own the status file (else we can only read it)
@var filename: Path for the status file (default to rc.statfile)
@var pid: Process ID for the status file owner
@var dataset: The dataset name (default to rc.dataset)
@var timestamp: The timestamp for the status file
@var progress: The number of steps taken
@var total: The total number of steps until finishing

@note: This is implemented as a module instead of an instantiable class,
because there is only supposed to be one status file active at a time.

                                   
"""

__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import logging as log
import os
import time

from mscanner.configuration import rc
import statusfile as s

owner = False
filename = None
pid = 0
dataset = ""
timestamp = 0
progress = 0
total = 0

def read(filename):
    """Attempt to read the status file. Returns True if successful, and False if
    the status file does not exist"""
    s.filename = filename if filename else rc.statfile
    if not s.filename.isfile():
        return False
    lines = s.filename.lines()
    s.owner = False
    s.pid = int(lines[0])
    s.progress = int(lines[1])
    s.total = int(lines[2])
    s.dataset = lines[3].strip()
    s.timestamp = float(lines[4])
    return True

def start(filename=None, total=0, dataset="", progress=0, timestamp=None):
    """Set the status file to beginning of operation.    
    """
    s.filename = filename if filename else rc.statusfile
    if s.filename.isfile() and not s.owner:
        raise IOError("Status file is in use!")
    s.pid = os.getpid()
    s.dataset = dataset if dataset else rc.dataset
    s.timestamp = timestamp if timestamp else time.time()
    s.owner = True
    s.update(progress, total, dolog=False)
    
def close():
    """Remove status file"""
    if s.owner:
        if s.filename is not None and s.filename.isfile():
            s.filename.remove()
        s.owner = False
        s.filename = None
        s.pid = 0
        s.dataset = ""
        s.timestamp = 0
        s.progress = 0
        s.total = 0

def contents():
    """Return status file contents"""
    return "%d\n%d\n%d\n%s\n%s\n" % (
        s.pid, s.progress, s.total, s.dataset, str(s.timestamp))

def update(progress=None, total=None, dolog=True):
    """Update progress of status file.  If progress is None, set to total."""
    if s.filename is not None and not s.owner:
        raise IOError("Cannot update status file which we do not own")
    if total is not None:
        s.total = total
    if progress is None:
        s.progress = s.total
    else:
        s.progress = progress
    if dolog:
        log.info("Completed %d out of %d", s.progress, s.total)
    if s.filename is not None:
        s.filename.write_text(contents())
