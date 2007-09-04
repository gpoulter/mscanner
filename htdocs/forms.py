"""High-level form construction, based on web.form by Aaron Swartz

@author: Aaron Swartz (modified by Graham Poulter)
"""

import copy
import re
import web
from web import utils, net

def attrget(obj, attr, value=None):
    if hasattr(obj, '__contains__') and attr in obj: return obj[attr]
    if hasattr(obj, attr): return getattr(obj, attr)
    return value

class Form:
    """
    @ivar inputs: List of input fields in the form
    @ivar valid: True if the form is unfilled, or validly filled
    @ivar note: Message about invalid stuff
    @ivar validators: List of validators that operate on the whole form
    """
    
    def __init__(self, *inputs, **kw):
        """Construct a form.  Positional parameters are the form inputs.
        
        @param validators: Optional keyword, providing, a list of additional
        validators on the form besides the ones associated with an input.  
        """
        self.inputs = inputs
        self.valid = True
        self.note = None
        self.validators = kw.pop('validators', [])


    def __call__(self, inputs=None):
        """To avoid overwriting, create a new form by calling the class"""
        newform = copy.deepcopy(self)
        if inputs: 
            newform.validates(inputs)
        return newform


    def render(self):
        """A bare rendering of the form"""
        out = ''
        if self.note: out += '<p class="error">'+self.note+'</p>\n'
        out += '<table class="form">\n'
        for i in self.inputs:
            out += '<tr class="input">\n'
            out += '<th>' + i.renderlabel() + '</th>\n'
            out += '<td class="value">' + i.pre+i.render()+i.post + '</td>\n'
            out += '</tr>\n'
            if i.note is not None:
                out += '<tr class="error">\n<td colspan="2">\n'
                out += i.note or ""
                out += '\n</td>\n</tr>\n'
        out += "</table>\n"
        return out
    
    
    def validates(self, source, _validate=True):
        """
        @param source: Storage object from which to retrieve values
        
        @returns: True/False about whether the form validates."""
        if hasattr(self, "d"): del self._d
        isvalid = True
        for i in self.inputs:
            value = attrget(source, i.name)
            if _validate:
                isvalid = isvalid and i.validate(value)
            else:
                i.value = value
        if _validate:
            isvalid = isvalid and self._validate(source)
            self.valid = isvalid
        return isvalid


    def _validate(self, value):
        """Run additional validators for the form"""
        self.value = value
        for v in self.validators:
            if not v.valid(value):
                self.note = v.msg
                return False
        return True


    def fill(self, source=None, **kw):
        """Fil the form with values"""
        return self.validates(source, _validate=False, **kw)
    
    
    def __getitem__(self, i):
        """Dictionary access to inputs"""
        for x in self.inputs:
            if x.name == i: return x
        raise KeyError, i
    
    
    @property
    def d(self):
        """A storage dictionary of the form inputs (deleted by validates())"""
        try:
            return self._d
        except AttributeError:
            self._d = utils.storage([(i.name, i.value) for i in self.inputs])
            return self._d
    
    
    @property
    def errors(self):
        """A storage dictionary of the form errors."""
        return utils.storage([(i.name, i.note) for i in self.inputs])



class Input(object):
    """Represents an input in a form
    
    Constructor:
    @ivar name: Name attribute for the input
    @ivar validators: List of validators for the input
    
    Constructor keywords:
    @ivar label: Contents of the <label>
    @ivar pre: Text before the input
    @ivar post: Text after the input
    @ivar id: For id= attribute (but defaults to name if not provided)
    @ivar attrs: Other attributes
    
    Derived:
    @ivar note: Set by the first validator that fails
    """
    
    def __init__(self, name, *validators, **attrs):
        """
        @note: Use keyword class_ to specify the class= attribute.
        """
        self.name = name
        self.note = None
        self.validators = validators
        self.label = attrs.pop('label', name)
        self.value = attrs.pop('value', None)
        self.pre = attrs.pop('pre', "")
        self.post = attrs.pop('post', "")
        self.id = attrs.setdefault("id", name)
        if 'class_' in attrs: 
            attrs['class'] = attrs['class_']
            del attrs['class_']
        self.attrs = attrs


    def validate(self, value):
        """Return true if all validators work, otherwise fals and sets note
        to the validator message"""
        self.value = value
        for v in self.validators:
            if not v.valid(value):
                self.note = v.msg
                return False
        return True


    def render(self):
        """Render the <input> element itself"""
        raise NotImplementedError
    
    
    def renderlabel(self):
        """Render the label for the input"""
        return '<label for="%s">%s</label>' % (self.id, self.label)


    def addatts(self):
        str = ""
        for (n, v) in self.attrs.items():
            str += ' %s="%s"' % (n, net.websafe(v))
        return str
    
    
    
class Textbox(Input):
    def render(self):
        x = '<input type="text" name="%s"' % net.websafe(self.name)
        if self.value: x += ' value="%s"' % net.websafe(self.value)
        x += self.addatts()
        x += '>'
        return x



class Password(Input):
    def render(self):
        x = '<input type="password" name="%s"' % net.websafe(self.name)
        if self.value: x += ' value="%s"' % net.websafe(self.value)
        x += self.addatts()
        x += '>'
        return x



class Textarea(Input):
    def render(self):
        x = '<textarea name="%s"' % net.websafe(self.name)
        x += self.addatts()
        x += '>'
        if self.value is not None: x += net.websafe(self.value)
        x += '</textarea>'
        return x



class Dropdown(Input):
    def __init__(self, name, args, *validators, **attrs):
        self.args = args
        super(Dropdown, self).__init__(name, *validators, **attrs)


    def render(self):
        x = '<select name="%s"%s>\n' % (net.websafe(self.name), self.addatts())
        for arg in self.args:
            if type(arg) == tuple:
                value, desc= arg
            else:
                value, desc = arg, arg 
            if self.value == value: 
                select_p = ' selected="selected"'
            else: 
                select_p = ''
            x += '  <option %s value="%s">%s</option>\n' % (
                select_p, net.websafe(value), net.websafe(desc))
        x += '</select>\n'
        return x



class Radio(Input):
    def __init__(self, name, args, *validators, **attrs):
        """
        @param args: (name, label) pairs for the radio buttons
        
        @note: The radio inputs do not get any attributes.
        """
        self.args = args
        super(Radio, self).__init__(name, *validators, **attrs)


    def renderlabel(self):
        """Dummy label: no unique ID for the set of buttons"""
        return self.label


    def render(self, only=None):
        """Write a list of radio inputs. If argname is a string, only write the
        radio button for that value."""
        out = ""
        for arg, label in self.args:
            if only is not None and arg != only:
                continue
            if label is None: label = arg
            if self.value == arg: 
                checked = ' checked="checked"'
            else: 
                checked = ''
            out += '<input type="radio" name="%s" value="%s"%s> %s ' % \
            (net.websafe(self.name), net.websafe(arg), checked, net.websafe(label))
        return out



class Checkbox(Input):
    def render(self):
        x = '<input name="%s" type="checkbox"' % net.websafe(self.name)
        if self.value: x += ' checked="checked"'
        x += self.addatts()
        x += '>'
        return x



class Button(Input):
    def __init__(self, name, *validators, **attrs):
        super(Button, self).__init__(name, *validators, **attrs)
        self.description = ""

    def render(self):
        safename = net.websafe(self.name)
        x = '<button name="%s"%s>%s</button>' % (safename, self.addatts(), safename)
        return x



class Hidden(Input):
    def __init__(self, name, *validators, **attrs):
        super(Hidden, self).__init__(name, *validators, **attrs)
        # it doesnt make sence for a hidden field to have description
        self.description = ""

    def render(self):
        x = '<input type="hidden" name="%s"' % net.websafe(self.name)
        if self.value: x += ' value="%s"' % net.websafe(self.value)
        x += ' />'
        return x



class File(Input):
    def render(self):
        x = '<input type="file" name="%s"' % net.websafe(self.name)
        x += self.addatts()
        x += '>'
        return x
    
    
    
class Validator:
    """Generic validator."""
    def __init__(self, test, msg, jstest=None): 
        """ 
        @param test: Function applied to input value when used as an input
        validator, and applied to the Storage source for the form when used as
        a form validator.
        
        @param msg: To be assigned to the note when validator fails.
        """
        utils.autoassign(self, locals())
        
    def __deepcopy__(self, memo): 
        return copy.copy(self)
    
    def valid(self, value): 
        try: 
            return self.test(value)
        except: 
            return False



class RegexValidator(Validator):
    """Tests that the value matches a particular regular expression"""
    def __init__(self, rexp, msg):
        self.rexp = re.compile(rexp)
        self.msg = msg
    
    def valid(self, value):
        return bool(self.rexp.match(value))


# Intended to match non-null inputs
notnull = Validator(
    bool, "Required")


# Matches if value is None or "ok"
checkbox_validator = Validator(
    lambda x: x == None or x == "on",
    "Bad checkbox")


def ischecked(value):
    """Boolean of a checkbox in web.input()"""
    if value == None:
        return False 
    elif value == "on":
        return True
    else:
        return bool(value)