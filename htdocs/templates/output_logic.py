import web
from . import output, query_logic
from htdocs import forms, helpers

OutputForm = forms.Form(
    forms.Hidden(
        "operation",
        forms.Validator(lambda x: x in ["download", "delete"], "Invalid op")),
    
    forms.Checkbox(
        "omit_mesh",
        forms.checkbox_validator),
    
    forms.Hidden(
        "dataset",
        forms.notnull,
        query_logic.dataset_validator),
    
    forms.Hidden(
        "delcode",
        query_logic.delcode_validator),
)



class OutputPage:
    """Page linking to outputs"""
    
    def GET(self):
        """Just list the available output directories"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        page = output.output()
        page.donetasks = helpers.list_visible()
        print page
        
        
    def POST(self):
        """Attempt to download or delete one of the outputs"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        page = output.output()
        page.donetasks = helpers.list_visible()
        oform = OutputForm()
        
        if not oform.validates(web.input()):
            e = ["<li>%s: %s</li>\n" % (n,e) for n,e in 
                 oform.errors.iteritems() if e is not None]
            page.errors = "".join(["<p>Errors</p><ul>\n"]+e+["</ul>\n"])
            page.inputs = oform
            print page
            return
            
        d = oform.d
        page.target = d.dataset
        if d.operation == "delete":
            status, descriptor = helpers.get_task_status(
                d.dataset, helpers.get_current_task())
            if status == "current":
                page.nodelete = "MScanner is busy with the task."
            elif status == "notfound":
                page.nodelete = "there is no task with that name."
            elif status in ["queue","done"]:
                if d.delcode != descriptor.delcode:
                    page.nodelete = "incorrect deletion code."
                else:
                    if status == "done":
                        # Attempt to remove the output directory
                        dirname = descriptor._filename.parent
                        try:
                            for fname in dirname.files():
                                fname.remove()
                            dirname.rmdir()
                        except (OSError, WindowsError), e:
                            # Whoopsie daisy - have permission?
                            print e
                        page.donetasks = helpers.list_visible()
                        page.deleted_location = "output"
                    elif target_status == "queue":
                        # Attempt to remove queue file for the target
                        descriptor._filename.remove()
                        page.deleted_location = "queue"
            print page
            
        elif d.operation == "download":
            # Save the output directory as a zip file for download
            from mscanner.configuration import rc
            dataset = d.dataset
            outdir = rc.web_report_dir / dataset
            outfile = outdir / (dataset + ".zip")
            if not outfile.exists():
                from zipfile import ZipFile, ZIP_DEFLATED
                zf = ZipFile(str(outfile), "w", ZIP_DEFLATED)
                for fpath in outdir.files():
                    if fpath.endswith(".zip"):
                        # Omit existing zip files
                       continue
                    if fpath.basename() == rc.report_term_scores:
                        # Omit MeSH terms if the user requests it
                        if forms.ischecked(oform.d.omit_mesh):
                            continue
                    if fpath.basename() == rc.report_result_all:
                        # Omit all-in-one result file
                        continue
                    zf.write(str(fpath), str(fpath.basename()))
                zf.close()
                outfile.chmod(0777)
            dsurl = web.urlquote(dataset)
            web.seeother("static/output/" + dsurl + "/" + dsurl + ".zip")
            print page
