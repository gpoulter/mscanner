#!/usr/bin/env python

"""Queueing facility for the web frontend

Queueing program checks queue directory every second for new descriptor files,
and starts a query or validation operation. When the operation completes, the
descriptor file is moved to the output.

Example descriptor file for query::
    #operation = query
    #dataset = Whatever
    #limit = 500
    #threshold = 10.3
    #submitted = 23424123.3
    804133
    3214241
    ...
    
Example descriptor file for validation::
    #operation = validate
    #dataset = Whatever
    #numnegs = 100000
    #alpha = 0.5
    #submitted = 23424123.3
    804133
    3214241
    ...
"""

from __future__ import with_statement
from __future__ import division

                                     
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

import logging
import os
from path import path
import sys
import time

from mscanner.configuration import rc
from mscanner.medline.Databases import Databases
from mscanner.core.QueryManager import QueryManager
from mscanner.core.ValidationManager import CrossValidation
from mscanner.core import iofuncs
from mscanner.scripts import update


def parsebool(s):
    """Handler for converting strings to booleans"""
    if isinstance(s, basestring):
        s = s.strip()
        if s == "0" or s == "False":
            return False
        elif s == "1" or s == "True":
            return True
        else:
            raise ValueError("Failed to parse boolean: %s" % s)
    else:
        return bool(s)


descriptor_keys = dict(
    alpha=float,
    code=str,
    dataset=str,
    delcode=str,
    hidden=parsebool,
    limit=int,
    numnegs=int,
    nfolds=int,
    operation=str,
    threshold=float,
    submitted=float,
    timestamp=float,)


def read_descriptor(fpath):
    """Reads a descriptor file, returning a dictionary of parameters.

    Each line is '#key = value'. We stops at the first line that not starting
    with '#'.  Valid keys are in L{descriptor_keys}. The same file can be used
    with read_pmids, which will ignores the lines beginning with '#'.

    @return: Storage object, with additional '_filename' containing fpath."""
    from mscanner.core.Storage import Storage
    result = Storage()
    with open(fpath, "r") as f:
        line = f.readline()
        while line.startswith("#"):
            key, value = line[1:].split(" = ",1)
            value = descriptor_keys[key](value.strip())
            result[key] = value
            line = f.readline()
        result["_filename"] = fpath
    return result


def write_descriptor(fpath, pmids, params):
    """Write parameters and PubMed IDs to the descriptor file.
    @param fpath: File to write
    @param pmids: List of PubMed IDs, may be None
    @param params: Dictionary to write. Values are converted with str(). Only
    keys from descriptor_keys are used."""
    with open(fpath, "w") as f:
        fpath.chmod(0777)
        for key, value in params.iteritems():
            if key in descriptor_keys: 
                f.write("#" + key + " = " + str(value) + "\n")
        if pmids is not None:
            for pmid in pmids:
                f.write(str(pmid)+"\n")


class QueueStatus:
    """Describes the current state of the queue
    
    @ivar tasklist: Descriptors of tasks in the queue, oldest first.
    
    @ivar running: First member of L{tasklist}, which is being processed.
    
    @ivar donelist: Completed tasks, oldest first.
    
    @ivar status: Mapping from dataset to status code (DONE, RUNNING, WAITING)
    
    @ivar _tasks: Mapping from dataset to task object
    """
    
    DONE = "done"
    RUNNING = "running"
    WAITING = "waiting"
    
    def __init__(self, with_done=True):
        """Constructor for the status
        
        @param with_done: Set this to False if you don't need L{donelist}."""
        self._load_tasklist()
        self.donelist = []
        if with_done: self._load_donelist()
        self._load_maps()

    
    def _load_tasklist(self):
        """Populate L{tasklist}.
        
        @note: We only load files that are older than one second.  Without
        this we sometimes catch files half-written by the web interface."""
        current_time = time.time()
        eligible_files = [f for f in rc.queue_path.files() \
                          if f.mtime < current_time-1]
        self.tasklist = [read_descriptor(f) for f in eligible_files]
        self.tasklist.sort(key=lambda x:x.submitted)
        self.running = self.tasklist[0] if self.tasklist else None


    def _load_donelist(self):
        """Populate L{donelist}"""
        self.donelist = []
        for fpath in rc.web_report_dir.dirs():
            if (fpath/rc.report_descriptor).exists():
                self.donelist.append(
                    read_descriptor(fpath/rc.report_descriptor))
        self.donelist.sort(key=lambda x:x.submitted)


    def _load_maps(self):
        """Calculate the L{status} and L{_tasks} mapping"""
        self.status = {}
        self._tasks = {}
        for task in self.tasklist:
            self.status[task.dataset] = self.WAITING
            self._tasks[task.dataset] = task
        for task in self.donelist:
            self.status[task.dataset] = self.DONE
            self._tasks[task.dataset] = task
        if self.tasklist:
            self.status[self.running.dataset] = self.RUNNING

    
    def __getitem__(self, dataset):
        """Retrieve the task descriptor for a given data set."""
        return self._tasks.__getitem__(dataset)


    def __contains__(self, dataset):
        """Return whether given dataset exists"""
        return self._tasks.__contains__(dataset)


    def position(self, dataset):
        """Return distance of dataset from front of queue."""
        for idx, d in enumerate(self.tasklist):
            if d.dataset == dataset:
                return idx
        return None

    
    

def delete_output(dataset):
    """Delete the output directory for the given task"""
    logging.debug("Attempting to delete output for %s" % dataset)
    dirpath = rc.web_report_dir / dataset
    for fname in dirpath.files():
        fname.remove()
    dirpath.rmdir()
    

def mainloop():
    """Look for descriptor files every second"""
    env = None
    try:
        last_clean = 0 # time.time() of last output-cleaning
        last_update = 0 # time.time() of last database update
        while True:
            # Cron: delete the oldest outputs
            if time.time() - last_clean > 6*3600:  
                logging.info("Looking for old datasets")
                queue = QueueStatus()
                queue.donelist.reverse() # Newest first
                for task in queue.donelist[100:]:
                    try:
                        delete_output(task.dataset)
                    except OSError:
                        pass # Failed to delete output
                last_clean = time.time()
            
            # Cron: update the databases
            if time.time() - last_update > 12*3600:
                if env is not None: env.close()
                env = None
                update.update_mscanner()
                env = Databases()
                env.article_list # long first load time
                last_update = time.time()
            
            # Now perform any queued tasks
            queue = QueueStatus()
            task = queue.running
            if task is not None:
                # The output directory for the task
                outdir = rc.web_report_dir / task.dataset
                logging.info("Starting %s for %s", task.operation, task.dataset)
                task._filename.utime(None) # Update mod time for status display
                try:
                    if task.operation == "query":
                        QM = QueryManager(
                            outdir, 
                            dataset=task.dataset,
                            limit=task.limit,
                            threshold=task.threshold, 
                            env=env)
                        QM.query(task._filename)
                        QM.write_report()
                        QM.__del__()
                    elif task.operation == "validate":
                        rc.alpha = task.alpha
                        VM = CrossValidation(
                            outdir, 
                            dataset=task.dataset,
                            env=env)
                        VM.validation(task._filename, task.numnegs)
                        VM.report_validation()
                        VM.__del__()
                    task._filename.move(outdir / "descriptor.txt")
                except ValueError, e:
                    logging.exception(e)
            else:
                # Wait before the next iteration
                time.sleep(1)
    finally:
        if env is not None: env.close()


def populate_test_queue():
    """Place some dummy queue files to test the queue operation"""
    from mscanner.core.Storage import Storage
    pmids = list(iofuncs.read_pmids(rc.corpora / "Test" / "gdsmall.txt"))
    task = Storage(
        operation = "validate", 
        dataset = "gdqtest_valid",
        numnegs = 1000, 
        alpha = 0.6,
        limit = 500, 
        threshold = 0, 
        submitted = time.time())
    write_descriptor(rc.queue_path/task.dataset, pmids, task)
    task.operation = "query"
    task.dataset = "gdqtest_query"
    task.submitted += 5
    write_descriptor(rc.queue_path/task.dataset, pmids, task)


if __name__ == "__main__":
    iofuncs.start_logger()
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        populate_test_queue()
    try:
        mainloop()
    except KeyboardInterrupt:
        pass
    logging.shutdown()
