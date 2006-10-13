#!/usr/bin/env python

"""PySlick Web Server

:Author: Graham Poulter
:Date: 2006/09/10
:Version: 20060910
:License: GPL-2
                                                             

Introduction
------------

PySlick is my idea of a pure-Python HTTPD capable of serving files and
processing .py CGI scripts.  It improves upon CGIHTTPServer by letting
you specify index files, an htdocs directory to serve out of, and by
working reliably through execfile().

Calls must be synchronous (cannot serve two requests while one is
still processing).  If you want to have multiple requests processing
at once (and have your CGI scripts know about locking, etc), use
Apache.

Usage
-----

Usage: ./pyslick.py <htdocs> <port>

<htdocs> is the webserver root (defaults to current directory). All
files under <htdocs> will be served.  <port> is the port to listen on
(defaults to 8080).

Functionality
-------------

 * Serves up files according to MIME type
 * Handles GET and POST requests
 * Supports configurable list of directory index files
 * Lists contents for directories without an index
 * Denies access to things outside self.htdocs
 * Uses execfile() on Python scripts, setting up CGI environment

Credits
-------

 * Started from the simple Python HTTPD by Jon Berg, found at
   turtlemeat.com

 * Used SimpleHTTPServer and CGIHTTPServer for ideas.

 * Modified CGIHTTPServer.CGIHTTPRequest.run_cgi to make the CGI
   environment in SlickRequestHandler.update_environment

Dependencies
------------

 * path module for Pythonic path manipulations.
   http://www.jorendorff.com/articles/python/path/

Contributions
-------------

 * If you want to add features or fix bugs, feel free to send me your
   modified code and/or a diff against this file, to graham dot
   poulter at gmail dot com.

"""

import cgi, mimetypes, os, sys, urllib
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer
from path import path
from cStringIO import StringIO

class SlickRequestHandler(BaseHTTPRequestHandler):
    """Request handler for the Slick HTTPD"""
    
    # Search list for the directory index
    indeces = ["index.py","index.xhtml","index.html","index.htm"]

    # Directory from which to serve documents
    htdocs = path(os.getcwd())

    # Mapping of extension to MIME type
    extensions_map = mimetypes.types_map.copy()
    extensions_map.update({
        '': 'application/octet-stream',
        '.py': 'text/plain',
        '.c': 'text/plain',
        '.h': 'text/plain',
        '.js': 'text/javascript',
        '.xhtml': 'text/xml',
        })

    def parse_path(self, path):
        """Process request URL, returning target and query strings

        Splits path into file and '?' query, the makes target relative
        to self.htdocs, denying access if result does not lie under
        it.  For directories, tries to obtain a directory index.
        """
        # Split around the "?"
        qmark = path.rfind("?")
        if qmark != -1:
            target, query = path[:qmark], path[qmark+1:]
        else:
            target, query = path, None
        # Target is relative to htdocs
        while target.startswith('/'):
            target = target[1:]
        target = (self.htdocs / target).normpath()
        # Target should be under self.htdocs
        if not target.startswith(self.htdocs):
                self.send_error(403, "Access denied to %s (not under %s)" % (target, self.htdocs))
                return None, None
        # Replace directory with index if possible
        if target.isdir():
            for index in self.indeces:
                if (target / index).isfile():
                    target = target / index
        return target, query

    def list_directory(self, target):
        """Helper to produce a directory listing (absent index.html).
        
        Prints the HTML with the directory listing (prints nothing if
        there is an error).
        """
        # Get contents of target directory
        os.chdir(self.htdocs)
        target = self.htdocs.relpathto(target)
        try:
            dirs = target.dirs()
            files = target.files()
        except os.error:
            self.send_error(404, "No permission to list directory")
            return
        dirs.sort(key=lambda a: a.lower())
        files.sort(key=lambda a: a.lower())
        # Write header
        f = StringIO()
        p1 = '-//W3C//DTD XHTML 1.0 Strict//EN'
        p2 = 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd'
        f.write('<!DOCTYPE html PUBLIC "%s" "%s">\n' % (p1,p2))
        f.write('<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n')
        f.write('<head>\n')
        f.write('<title>Directory listing for %s</title>\n' % target)
        f.write('</head>\n')
        f.write('<body>\n')
        f.write('<h2>Directory listing for %s</h2>\n' % target)
        f.write('<hr/>\n')
        # Write listing
        f.write('<a href="/%s">[PARENT]</a>\n' % urllib.quote(target.parent))
        f.write('<ul>\n')
        for dirname in dirs:
            f.write('<li><a href="/%s">%s</a></li>\n' % (urllib.quote(dirname), cgi.escape(dirname.basename()+"/")))
        for fname in files:
            displayname = fname.basename()
            if fname.islink(): displayname += "@"
            f.write('<li><a href="/%s">%s</a></li>\n' % (urllib.quote(fname), cgi.escape(displayname)))
        f.write('</ul>\n')
        f.write('<hr/>\n')
        f.write('</body>\n')
        f.write('</html>\n')
        length = f.tell()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(length))
        self.end_headers()
        self.wfile.write(f.getvalue())

    def update_environment(self, target, query=None):
        """Set environment variables for the CGI script.

        Update os.environ with server, path, query, host, remote,
        authorization, content-type, content-length, user-agent and
        cookie information.
        """
        env = dict()
        # Server info
        env["SERVER_SOFTWARE"] = self.version_string()
        env["SERVER_NAME"] = self.server.server_name
        env["GATEWAY_INTERFACE"] = "CGI/1.1"
        env["SERVER_PROTOCOL"] = self.protocol_version
        env["SERVER_PORT"] = str(self.server.server_port)
        env["REQUEST_METHOD"] = self.command
        env["SCRIPT_NAME"] = target.relpath()
        # Path to script
        env["PATH_INFO"] = urllib.unquote(target.relpath())
        # Query string
        if query is not None:
            env["QUERY_STRING"] = query
        # Host address
        host = self.address_string()
        if host != self.client_address[0]:
            env["REMOTE_HOST"] = host
        # Remote address
        env["REMOTE_ADDR"] = self.client_address[0]
        # Authorization
        authorization = self.headers.getheader("Authorization")
        if authorization:
            authorization = authorization.split()
            if len(authorization) == 2:
                import base64, binascii
                env["AUTH_TYPE"] = authorization[0]
                if authorization[0].lower() == "basic":
                    try:
                        authorization = base64.decodestring(authorization[1])
                    except binascii.Error:
                        pass
                    else:
                        authorization = authorization.split(":")
                        if len(authorization) == 2:
                            env["REMOTE_USER"] = authorization[0]
        # content-type header
        if self.headers.typeheader is None:
            env["CONTENT_TYPE"] = self.headers.type
        else:
            env["CONTENT_TYPE"] = self.headers.typeheader
        # content-length header
        length = self.headers.getheader("Content-Length")
        if length:
            env["CONTENT_LENGTH"] = length
        # accept header
        accept = []
        for line in self.headers.getallmatchingheaders("Accept"):
            if line[:1] in "\t\n\r ":
                accept.append(line.strip())
            else:
                accept = accept + line[7:].split(",")
        env["HTTP_ACCEPT"] = ",".join(accept)
        # user-agent header
        ua = self.headers.getheader("User-Agent")
        if ua:
            env["HTTP_USER_AGENT"] = ua
        # cookie header
        co = filter(None, self.headers.getheaders("Cookie"))
        if co:
            env["HTTP_COOKIE"] = ", ".join(co)
        # provide defaults for various headers
        for k in ("QUERY_STRING", "REMOTE_HOST", "CONTENT_LENGTH", "HTTP_USER_AGENT", "HTTP_COOKIE"):
            env.setdefault(k, "")
        # update environment
        os.environ.update(env)

    def run_cgi(self, target, query=None):
        """Run a Python CGI script.

        First sets the environment, saves stdout etc., does
        execfile(), and restores stdout etc.
        """
        # functions to save/restore argv/stdin/stdout/stderr
        def get_sys():
            return sys.argv, sys.stdin, sys.stdout, sys.stderr
        def set_sys(args):
            sys.argv, sys.stdin, sys.stdout, sys.stderr = args
        if self.htdocs not in sys.path:
            sys.path.insert(0, self.htdocs)
        self.update_environment(target, query)
        self.send_response(200, "Script output follows")
        print "HEADERS: " + str(self.headers)
        #print "INPUT: " + self.rfile.read(int(self.headers["content-length"]))
        #return
        try:
            saves = get_sys()
            os.chdir(target.dirname())
            try:
                set_sys(([target], self.rfile, self.wfile, self.wfile))
                execfile(target, {"__name__": "__main__"})
            finally:
                set_sys(saves)
        except SystemExit, sts:
            self.log_error("CGI script exit status %s", str(sts))
        else:
            self.log_message("CGI script exited OK")

    def do_GET(self):
        """Process an HTTP GET request

        Gives back file contents based on MIME type.  For directories
        without an index, lists contents.  For .py files, executes as
        CGI scripts.
        """
        target, query = self.parse_path(self.path)
        if target is None:
            return
        # Directory listing
        if target.isdir():
            self.list_directory(target)
        # CGI script for .py files
        elif target.ext == ".py":
            self.log_message("Target is %s" % target)
            self.run_cgi(target, query)
        # Return file according to MIME type
        else:
            if not target.isfile():
                self.send_error(404, "File not found: %s" % target)
                return
            self.send_response(200)
            if target.ext in self.extensions_map:
                self.send_header('Content-Type', self.extensions_map[target.ext.lower()] + '; charset=utf-8')
            else:
                self.send_header('Content-Type', self.extensions_map[''] + '; charset=utf-8;')
            self.end_headers()
            try:
                self.wfile.write(target.bytes())
            except IOError:
                self.send_error(404, "File unreadable: %s" % target)

    def do_POST(self):
        """Process an HTTP POST request.

        The request must be to a .py file, which is called as a CGI
        script"""
        target, query = self.parse_path(self.path)
        if target is None:
            return
        if not target.isfile():
            self.send_error(404, "File not found: %s" % target)
            return
        if target.ext != ".py":
            self.send_error(403, "Unrecognised extension: %s" % target.ext)
            return
        self.run_cgi(target, query)
        
def serverLoop(htdocs=None, port=8080):
    """Run a PySlick HTTPD on localhost at the given port"""
    if htdocs is not None:
        SlickRequestHandler.htdocs = htdocs
    try:
        server = HTTPServer(('',port), SlickRequestHandler)
        print 'Starting PySlick HTTPD...'
        server.serve_forever()
    except KeyboardInterrupt:
        print '^C received, shutting down PySlick HTTPD'
        server.socket.close()

if __name__ == '__main__':
    htdocs = None
    port = 8080
    if len(sys.argv) >= 2:
        htdocs = path(sys.argv[1]).abspath()
        if not htdocs.isdir():
            raise ValueError("htdocs path %s does not exist" % htdocs)
    if len(sys.argv) >= 3:
        port = int(sys.argv[2])
    serverLoop(htdocs, port)
