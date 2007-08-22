#!/export/home/medscan/local32/bin/python2.5

"""
web.py controller for the MScanner web interface (view is the
template code, and the model is the queue programme).

                                   
"""

__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""

import sys
sys.path.insert(0,"/export/home/medscan/source")

import web
from web.utils import Storage
import re
import time

from templates import contact, form, front, form, output, status
from mscanner.configuration import rc
from mscanner.scorefile import readDescriptor, writeDescriptor

# Set pretty and informative error handler
web.webapi.internalerror = web.debugerror

# Configure URLs for the application
urls = (
    '/', 'web_front',
    '/form', 'web_form',
    '/status', 'web_status',
    '/output', 'web_output',
    '/contact', 'web_contact',
)

# Some variables that all pages use
page_defaults = dict(
    baseurl = "",
    siteurl = rc.siteurl,
    inc = lambda x: "templates/" + x,
    )

# Default values to fill in on the form
form_defaults = Storage(
    alpha = 0.5,
    captcha = "",
    delcode = "",
    dataset = "",
    limit = 500,
    numnegs = 50000,
    nfolds = 10,
    operation = "",
    positives = "",
    threshold = 20.0,
    )

#############################################################
## Helper functions
#############################################################

def makeQueueFilename(timestamp, dataset):
    """Given time and task name, return the path for the queue file"""
    ftime = time.strftime("%Y%m%d-%H%M%S", time.gmtime(timestamp))
    return rc.queue_path / ("%s-%s.txt" % (ftime, dataset))

def getStatusForTask(dataset, current=None):
    """Figure out the status of the given task identifier"""
    if current is not None and dataset == current.dataset:
        # Task is the one being processed
        return "current", current
    elif (rc.web_report_dir / dataset / rc.report_index).exists():
        # Task is is completed, descriptor is in the output directory
        return "done", readDescriptor(
            rc.web_report_dir / dataset / rc.report_descriptor)
    else:
        for queue_dataset, fname in zip(*listQueueTasks()):
            if dataset == queue_dataset:
                # Task is in the queue
                return "queue", readDescriptor(rc.queue_path / fname)
        # Task is nowhere to be found
        return "notfound", None

def checkDataset(dataset):
    """For security, we always check the task name is valid"""
    if re.match(r"^[a-zA-Z0-9._-]+$", dataset) is None:
        raise ValueError("Invalid task name %s!" % dataset)

def listQueueTasks(descriptors=False):
    """Return a list of tasks currently in the queue and a corresponding
    list of the filenames.

    @param descriptors: If True, return full descriptor objects instead
    of just the task name. """
    tasks = []
    filenames = rc.queue_path.files("*.txt")
    for filename in filenames:
        fmatch = re.match(r"\d{8}-\d{6}-(.+)\.txt", filename.basename())
        if fmatch is not None:
            if descriptors == True:
                tasks.append(readDescriptor(filename))
            else:
                tasks.append(fmatch.group(1))
    return tasks, filenames

def listDoneTasks(descriptors=False):
    """Return a list of completed task names. 
    
    @param descriptors: If True, return full descriptor objects instead
    of just the task name. """
    output = []
    for fpath in rc.web_report_dir.dirs():
        if (fpath / rc.report_index).exists():
            if descriptors == True:
                output.append(readDescriptor(fpath / rc.report_descriptor))
            else:
                output.append(fpath.basename())
    return output

def getCurrentTask():
    """Returns None if the queue is empty.  If busy, returns a Storage
    descriptor for the currently running task, with added queue_length
    variable for total number of tasks in the queue."""
    listing = rc.queue_path.files("*.txt")
    if not listing: return None
    current = readDescriptor(listing[0])
    current.queue_length = len(listing)
    return current
    
def validateForm(input):
    """Given a Storage object with keys as in form_defaults, return a
    dictionary of converted values, and a list of error strings to
    display if conversion was unsucessful on some of them"""
    validator = dict(
        alpha=(
            float, lambda x: 0.0 <= x <= 1.0, 
            "Alpha must be between 0.0 and 1.0"),
        captcha=(
            str, lambda x: x == "orange", 
            "The captcha must be the word 'orange'"),
        delcode=(
            str, lambda x: re.match(r"^[ a-zA-Z0-9.;_-]{0,20}$", x) is not None,
            "Deletion code must be alphanumeric and shorter than 20 characters"),
        dataset=(
            lambda x: re.sub(r"[^ a-zA-Z0-9.;_-]", "", x)[:50], 
            lambda x: re.match(r"^[a-zA-Z0-9._-]{1,50}$", x) is not None,
            "The task name must contain a-z, ., _, - characters only "+
            "as it needs to be safe as a directory name in a URL."),
        limit=(
            int, lambda x: 100 <= x <= 10000,
            "Retrieval limit must be between 100 and 10000."),
        nfolds=(
            int, lambda x: x == 5 or x == 10,
            "Number of cross validation folds must be 5 or 10."),
        numnegs=(
            int, lambda x: 100 <= x <= 100000,
            "Number of negatives must be between 100 and 100000."),
        operation=(
            str, lambda x: x in ["delete", "download", "query",  "validate"],
            "Operation must be 'delete', 'download', 'query' or 'validate'"),
        positives=(
            lambda x: [int(y) for y in x.split()], 
            lambda x: len(x) > 0, 
            "At least one PubMed ID must be provided"),
        threshold=(
            float, lambda x: -1000 <= x <= 1000,
            "The threshold must be between -1000 and +1000"),
        timestamp=(
            float, lambda x: x > 0,
            "Time stamp must be a floating point number.")
        )
    # List of errors
    errors = []
    # Place for output variables
    output = Storage()
    for key, value in input.iteritems():
        if key not in validator:
            output[key] = value
            continue
        convertfn, test, errmsg = validator[key]
        # Don't want to write all of the data if it's long
        svalue = str(value)
        errdata = svalue if len(svalue) < 25 else svalue[:25] + "..."
        # Error includes variable name, value, and error msg
        errmsg = "Error in " + key + " (" + errdata + "): " + errmsg
        if value is None:
            output[key] = None
        else:
            try:
                output[key] = convertfn(value) 
                if not test(output[key]):
                    errors.append(errmsg)
            except:
                errors.append(errmsg)
    return output, errors

#############################################################
## Classes that handle requests
#############################################################

class web_front:
    """Front page of the site"""
    
    def GET(self):
        """Return the front page for MScanner"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        t = front.front(searchList=[page_defaults])
        t.current = getCurrentTask()
        print t

class web_form:
    """Main submission form for queries or validation"""
    
    def GET(self):
        """Return the form with default values"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        t = form.form(searchList=[page_defaults, form_defaults])
        t.current = getCurrentTask()
        print t
        
    def POST(self):
        """Submit the form.  If errors, return it with the same values
        as it was submitted with, but with errors listed at the top.  
        Without errors, place an entry in the queue and go to the
        special status page for that task.
        """
        web.header('Content-Type', 'text/html; charset=utf-8') 
        import pprint
        pp = pprint.PrettyPrinter()
        # Submitted form must have EVERY key or you get "bad request"
        input = web.input(*form_defaults.keys())
        # Check prior validity of the input
        formvars, errors = validateForm(input)
        # Check that data set does not already exist
        if formvars.dataset in listQueueTasks()[0]:
            errors.append("The specified task name is already in the queue")
        elif formvars.dataset in listDoneTasks():
            errors.append("The specified task name is already in the output directory")
        #print pp.pformat(formvars)
        #print pp.pformat(errors)
        #return
        if errors:
            # Failed to convert everything - show the submission form again
            t = form.form(searchList=[page_defaults, input])
            t.current = getCurrentTask()
            t.errors = errors
            print t
        else:
            # All success - create a descriptor in the queue
            formvars.timestamp = time.time()
            destfile = makeQueueFilename(formvars.timestamp, formvars.dataset)
            writeDescriptor(destfile, formvars.positives, formvars)
            # Now go see the status page for this task
            web.seeother("status?dataset=%s;delcode=%s" % 
                         (formvars.dataset, web.urlquote(formvars.delcode)))

class web_status:
    """Status page.  
    
    Contains the current operation (if any), the status of the
    operation specified in the 'dataset' and 'delcode' variables,
    the list of tasks in the queue, and the contents of the log file.
    """
    
    def GET(self):
        """Output the status page"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        # Fill the general status template
        t = status.status(searchList=[page_defaults])
        t.current = getCurrentTask()
        t.listing = [f.basename() for f in rc.queue_path.files("*.txt")]
        t.logcontents = rc.logfile.lines()[-30:]
        # Check if we've been provided a data set to check on
        input = web.input(dataset=None, delcode=None)
        formvars, errors = validateForm(input)
        if errors:
            # Invalid dataset or delcode
            t.errors = errors
        elif formvars.dataset:
            # Give status for a particular dataset
            checkDataset(formvars.dataset)
            t.target = formvars.dataset 
            t.target_status, t.target_descriptor = \
             getStatusForTask(formvars.dataset, t.current)
            # Provide any deletion code to the form
            if formvars.delcode:
                t.delcode = formvars.delcode
        print t

class web_output:
    """Page linking to outputs"""
    
    def GET(self):
        """Just list the available output directories"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        t = output.output(searchList=page_defaults)
        t.donetasks = sorted(listDoneTasks(descriptors=True), 
                             key=lambda x:x.timestamp, reverse=True)
        print t
        
    def POST(self):
        """Attempt to download or delete one of the outputs"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        t = output.output(searchList=page_defaults)
        t.donetasks = listDoneTasks(descriptors=True)
        input = web.input("operation", dataset=None, delcode=None, omit_mesh=[])
        formvars, errors = validateForm(input)
        # So the returned form knows what task is referred to
        if formvars.dataset is None:
            errors.append("No data set was selected for the operation")
        else:
            t.target = formvars.dataset
        if errors:
            # There were invalid fields, so never mind operation, just print
            t.errors = errors
        elif formvars.operation == "delete":
            # Attempt to delete a specified task
            if formvars.delcode is None:
                # Must have deletion code to do deletion
                t.nodelete = "nodelcode"
            # Get status of the task we're supposed to delete
            target_status, target_descriptor = getStatusForTask(
                formvars.dataset, getCurrentTask())
            if target_status == "current":
                # Can't delete tasks in-progress
                t.nodelete = "busy"
            elif target_status == "notfound":
                # Can't delete tasks that don't exist
                t.nodelete = "notfound"
            elif target_status in ["queue","done"]:
                # Can delete if it's in the queue or finished
                if formvars.delcode != target_descriptor.delcode:
                    # Can't delete if the code's don't match
                    t.nodelete = "badcode"
                else:
                    if target_status == "done":
                        # Attempt to remove the output directory
                        dirname = target_descriptor._filename.parent
                        try:
                            for fname in dirname.files():
                                fname.remove()
                            dirname.rmdir()
                        except (OSError, WindowsError), e:
                            # Will deal with these better later
                            print e
                        t.donetasks = listDoneTasks(descriptors=True)
                        t.did_delete = "output"
                    elif target_status == "queue":
                        # Attempt to remove queue file for the target
                        target_descriptor._filename.remove()
                        t.did_delete = "queue"
        elif formvars.operation == "download":
            # Save the output directory as a zip file for download
            dataset = formvars.dataset
            outdir = rc.web_report_dir / dataset
            outfile = outdir / (dataset + ".zip")
            if not outfile.exists():
                from zipfile import ZipFile, ZIP_DEFLATED
                zf = ZipFile(str(outfile), "w", ZIP_DEFLATED)
                for fpath in outdir.files():
                    if fpath.endswith(".zip"):
                       continue
                    if fpath.basename() == rc.report_term_scores:
                        # Omit MeSH terms if the user requests it
                        if formvars.omit_mesh == ["on"]:
                            continue
                    zf.write(str(fpath), str(fpath.basename()))
                zf.close()
            dsurl = web.urlquote(dataset)
            web.seeother("static/output/" + dsurl + "/" + dsurl + ".zip")
        if formvars.operation != "download" or hasattr(t, "errors"):
            print t

class web_contact:
    """Form to contact the webmaster"""
    
    def GET(self):
        """Print the contact form"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        t = contact.contact(searchList=[page_defaults])
        print t
    
    def POST(self):
        """Submit the contact form.  Errors out if the captcha is 
        incorrect"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        ivars = web.input("captcha", "email", "message", "name")
        if ivars.email == "":
            ivars.email = "nobody@maples.stanford.edu"
        v = Storage()
        v.update(ivars)
        t = contact.contact(searchList=[page_defaults, v])
        if ivars.captcha.lower() != "orange":
            # Incorrect captcha
            v.bad_captcha = True
        else:
            # Send the email using the ivars
            import smtplib
            from email.mime.text import MIMEText
            msg = MIMEText(ivars.message)
            msg['Subject'] = "Contact from MScanner"
            msg['From'] = ivars.name + "<" + ivars.email + ">"
            msg['To'] = rc.webmaster_email
            server = smtplib.SMTP(rc.smtpserver)
            try:
                server.sendmail(ivars.email, rc.webmaster_email, msg.as_string())
            except Exception, e:
                v.failed_send = True
                v.error = str(e)
                sys.stderr.write(e.message)
            else:
                v.success = True
            server.quit()
        print t

if __name__ == "__main__":
    try:
        web.run(urls, globals())
    except KeyboardInterrupt:
        pass
