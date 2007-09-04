"""Logic behind the submission form
"""

import web

from . import query
from htdocs import helpers, forms
from mscanner.configuration import rc
from mscanner import scorefile

def parse_pmids(pmids):
    return [int(y) for y in pmids.split()]

delcode_validator = forms.RegexValidator(
    r"^[ a-zA-Z0-9.;:_-]{0,10}$", 
    "Should be 0-10 characters long, containing "+
    "only letters, numbers and .;:,_- punctuation.")

dataset_validator = forms.RegexValidator(
    r"^[ a-zA-Z0-9.,;:_-]{1,30}$",
    "Should be 1-30 characters long, containing "+
    "only letters, numbers and .,;:_- punctuation.")


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
        forms.Validator(
            lambda x: len(parse_pmids(x)) > 0,
            "Should be numbers separated by line breaks"),
        label="Input Citations", rows=3, cols=10),
          
    forms.Textbox(
        "dataset", dataset_validator, 
        forms.Validator(
            lambda x: x not in helpers.list_queue()[0], "Already in the queue"),
        forms.Validator(
            lambda x: x not in helpers.list_done(), "Already in the output"),
        label="Task Name", size=30),
          
    forms.Textbox(
        "delcode", delcode_validator, label="Deletion Code", size=8),
          
    forms.Checkbox(
        "hidden", forms.checkbox_validator, label="Hide output"),
    
    forms.Textbox(
        "threshold", 
        forms.Validator(
            lambda x: -1000 <= float(x) <= 1000,
            "Should be between -1000 and +1000"),
        label="Score Threshold", size=8),
    
    forms.Textbox(
        "limit", 
        forms.Validator(
            lambda x: 100 <= int(x) <= 10000,
            "Should be between 100 and 10000."),
        label="Result Limit", size=8),
          
    forms.Textbox(
        "numnegs", 
        forms.Validator(
            lambda x: 100 <= int(x) <= 100000,
            "Should be between 100 and 100000."),
        label="Number of Negatives", size=8),
          
    forms.Textbox(
        "alpha", 
        forms.Validator(
            lambda x: 0.0 <= float(x) <= 1.0,
            "Should be between 0.0 and 1.0"),
        label="F-Measure Alpha", size=8),
          
    forms.Radio(
        "operation",
        [ ("query", "Query Operation"), 
          ("validate", "Validation Operation") ],
        forms.Validator(lambda x: x in ["query", "validate"], "Invalid operation")),
)


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


class QueryPage:
    """Submission form for queries or validation"""
    
    def GET(self):
        """Return the form with default values"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        qform = QueryForm()
        qform.fill(form_defaults)
        page = query.query()
        page.inputs = qform
        page.current = helpers.get_current_task()
        print page
        
    def POST(self):
        """Submit the query form."""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        qform = QueryForm()
        if qform.validates():
            # Add a descriptor to the queue
            d = qform.d
            import time; d.timestamp = time.time()
            d.hidden = forms.ischecked(d.hidden)
            fpath = helpers.make_queue_name(d.timestamp, d.dataset)
            scorefile.writeDescriptor(fpath, parse_pmids(d.positives), d)
            # Status page for the task
            web.seeother("status?dataset=%s;delcode=%s" % 
                         (d.dataset, web.urlquote(d.delcode)))
        else:
            # Errors in the form, print it again
            page = query.query()
            page.inputs = qform
            page.current = helpers.get_current_task()
            print page
