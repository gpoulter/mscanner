"""web.py handler for the contact page"""

import web
import contact
from mscanner.htdocs import forms

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

ContactForm = forms.Form(
    
    #forms.Textbox(
    forms.Hidden(
        "captcha",
        forms.Validator(lambda x: x == "orange", "Should be the word 'orange'"),
        label="The word 'orange'", size=10),
    
    forms.Textbox(
        "name",
        forms.Validator(lambda x: len(x) < 50, "Should be less than 50 characters"),
        label="Name (optional)", size=35),
    
    forms.Textbox(
        "email",
        forms.Validator(lambda x: len(x) < 80, "Should be less than 80 characters"),    
        label="Email (optional)", size=35), 
    
    forms.Textarea(
        "message",
        forms.Validator(lambda x: len(x) < 2000, "Should be less than 2000 characters"),
        label="Message", rows=10, cols=40),

)
"""Structure for the form on the contact page"""



class ContactPage:
    """Form to contact the webmaster"""
    
    def GET(self):
        """Print the contact form"""
        web.header('Content-Type', 'text/html; charset=utf-8') 
        page = contact.contact()
        page.inputs = ContactForm()
        page.inputs["captcha"].value = "orange"
        print page
    
    
    def POST(self):
        """Submit the contact form."""
        import re
        sane = lambda s: re.sub(r"[/!#$%^&*{}[]|\\]+", "", s)
        web.header('Content-Type', 'text/html; charset=utf-8') 
        page = contact.contact()
        cform = ContactForm()
        if cform.validates(web.input()):
            try:
                import smtplib
                from email.mime.text import MIMEText
                from mscanner.configuration import rc
                fromaddr = sane(cform.d.email) or "nobody@maples.stanford.edu"
                msg = MIMEText(sane(cform.d.message))
                msg['Subject'] = "Contact from MScanner"
                msg['From'] = sane(cform.d.name) + " <" + fromaddr + ">"
                msg['To'] = rc.webmaster_email
                server = smtplib.SMTP(rc.smtpserver)
                server.sendmail(fromaddr, rc.webmaster_email, msg.as_string())
            except Exception, e:
                page.success = False
                page.error = str(e)
                page.inputs = cform
            else:
                page.success = True
                page.inputs = ContactForm()
            server.quit()
            page.inputs["captcha"].value = "orange"
            print page
        else:
            page.inputs = cform
            page.inputs["captcha"].value = "orange"
            print page
