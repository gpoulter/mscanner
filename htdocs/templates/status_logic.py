"""web.py handler for the status page"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import web
import status, query_logic
from htdocs import helpers, forms
from mscanner.configuration import rc


StatusForm = forms.Form(
    forms.Textbox(
        "dataset",
        query_logic.dataset_validator,
        label="Task name"
        ),
    forms.Textbox(
        "delcode",
        query_logic.delcode_validator,
        label="Deletion code"
        ),
)
"""Structure for the delete-this-task form on the status page"""


class StatusPage:
    """Lists the current status of MScanner and a given task.
    
    If the dataset and delcode parameters are given over the web,
    it provides a form for deleting the specified task.
    """
    
    def GET(self):
        """Print the status page"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        page = status.status()
        # Descriptor of currently running task (or None)
        page.current = helpers.get_current_task()
        # List of files in the queue
        page.listing = [f.basename() for f in rc.queue_path.files("*.txt")]
        # Last 30 lines of the log file
        page.logcontents = rc.logfile.lines()[-30:]
        inputs = web.input()
        if "dataset" in inputs:
            # Offer a form for deleting outputs
            if "delcode" not in inputs:
                inputs.delcode = ""
            sform = StatusForm(inputs)
            if sform.valid:
                page.inputs = sform
                page.target = sform.d.dataset
                page.target_status, page.target_descriptor = \
                 helpers.get_task_status(page.target, page.current)
        print page
