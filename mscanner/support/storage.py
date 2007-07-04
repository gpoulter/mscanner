"""Dictionary subclasses supporting dotted access

Storage -- Dictionary which supports d.foo attribute access
RCStorage -- Storage where d.foo returns d.foo() if d['foo'] is callable

                                   
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

class Storage(dict):
    """A Storage object is a dictionary, except `obj.foo` can be used in
    addition to `obj['foo']`.  Raises AttributeError instead
    of KeyError when dotted access is used."""
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError, k:
            raise AttributeError, k

    def __setattr__(self, key, value): 
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError, k:
            raise AttributeError, k

    def __str__(self):
        return "Storage(\n" + \
               "\n".join("   " + k + " = " + repr(v) + "," 
                         for k, v in self.iteritems()) + "\n)"

    def __repr__(self):     
        return '<Storage ' + dict.__repr__(self) + '>'

class RCStorage(Storage):
    """Subclass of Storage which auto-calls any callables 
    retrieved by dotted syntax (rc.foo is equivalent to
    rc.foo() if rc.foo is callable)

    This is intended for 'rc.foo = lambda: rc.bar + 10' statements
    so the value of rc.foo depends on the current value of rc.bar.
    """
    def __getattr__(self, key):
        v = Storage.__getattr__(self, key)
        return (v() if hasattr(v, "__call__") else v)

    def __str__(self):
        result = "RCStorage("
        for k, v in self.iteritems():
            if hasattr(v, "__call__"):
                v = v()
            result += "\n   " + k + " = " + str(v) + ","
        return result+"\n)"
