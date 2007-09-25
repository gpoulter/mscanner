"""web.py handler for the query submission page"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

import web
import time
import md5

import query
from htdocs import forms, queue
from mscanner.configuration import rc


def parse_pmids(pmids):
    """Parse a string into a list of integer PubMed IDs"""
    return [int(y) for y in pmids.split()]


delcode_validator = forms.RegexValidator(
    r"^[ a-zA-Z0-9.;:_-]{0,10}$", 
    "Should be 0-10 characters long, containing "+
    "only letters, numbers and .;:,_- punctuation.")
"""Checks deletion code for valid format"""


dataset_validator = forms.RegexValidator(
    r"^[ a-zA-Z0-9.,;:_-]{1,30}$",
    "Should be 1-30 characters long, containing "+
    "only letters, numbers and .,;:_- punctuation.")
"""Checks task name for valid format"""


def task_does_not_exist(dataset):
    """True if task does not exist in queue or output directory"""
    return not task_exists(dataset)

def task_exists(dataset):
    """True if task exists in queue or output directory"""
    return (rc.queue_path / dataset).isfile() or\
           (rc.web_report_dir / dataset).isdir()


QueryForm = forms.Form(
    forms.Hidden(
        "captcha",
        forms.Validator(
            lambda x: x == "orange", "Should be the word 'orange'"),
        label="Enter the word 'orange'"),
    
    forms.Hidden(
        "nfolds", 
        forms.Validator(
            lambda x: int(x) in [5,10], "Should be 5 or 10."),
        label="Number of validation folds"),
    
    forms.Textarea(
        "positives", 
        forms.Validator(lambda x: len(parse_pmids(x)) > 0,
            "Should be numbers separated by line breaks"),
        label="Input Citations", rows=3, cols=10),
    
    forms.Textbox(
        "dataset", 
        dataset_validator, 
        forms.Validator(task_does_not_exist, "Task already exists"),
        label="Task Name", size=30),
    
    forms.Textbox(
        "delcode", delcode_validator, label="Deletion Code", size=8),
    
    forms.Checkbox(
        "hidden", forms.checkbox_validator, label="Hide output"),
    
    forms.Textbox(
        "threshold", 
        forms.Validator(lambda x: -100 <= float(x) <= 100,
            "Should be between -1000 and +1000"),
        label="Score Threshold", size=8),
    
    forms.Textbox(
        "limit", 
        forms.Validator(lambda x: 100 <= int(x) <= 10000,
            "Should be between 100 and 10000."),
        label="Result Limit", size=8),
    
    forms.Textbox(
        "numnegs", 
        forms.Validator(lambda x: 100 <= int(x) <= 100000,
            "Should be between 100 and 100000."),
        label="Number of Negatives", size=8),
    
    forms.Textbox(
        "alpha", 
        forms.Validator(lambda x: 0.0 <= float(x) <= 1.0,
            "Should be between 0.0 and 1.0"),
        label="F-Measure Alpha", size=8),
    
    forms.Radio(
        "operation",
        [ ("query", "Query Operation"), 
          ("validate", "Validation Operation") ],
        forms.Validator(lambda x: x in ["query", "validate"], 
                        "Invalid operation")),
)
"""Structure of the query form"""


# Initial values to fill into the form
form_defaults = dict(
    alpha = 0.5,
    captcha = "orange",
    delcode = "",
    dataset = "",
    limit = 500,
    hidden = False,
    numnegs = 50000,
    nfolds = 10,
    operation = "query",
    positives = "",
    threshold = 20.0,
)
"""Default values for the query form"""



class QueryPage:
    """Submission form for queries or validation"""
    
    def GET(self):
        """Print the query form, filled with default values"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        page = query.query()
        page.inputs = QueryForm()
        page.inputs.fill(form_defaults)
        print page


    def POST(self):
        """Submit the query form, maintains previous values"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        qform = QueryForm()
        if qform.validates(web.input()):
            # Add a descriptor to the queue
            inputs = qform.d
            inputs.submitted = time.time()
            inputs.hidden = forms.ischecked(inputs.hidden)
            delcode_plain = inputs.delcode
            inputs.delcode = md5.new(delcode_plain).hexdigest()
            queue.write_descriptor(rc.queue_path / inputs.dataset, 
                                   parse_pmids(inputs.positives), inputs)
            # Show status page for the task
            web.seeother("status?dataset=%s;delcode=%s" % 
                         (inputs.dataset, web.urlquote(delcode_plain)))
        else:
            # Errors in the form, print it again
            page = query.query()
            page.queue = queue.QueueStatus(with_done=True)
            page.inputs = qform
            print page
