#!/usr/bin/python
import sys

print("sys.version[0]: %s" % sys.version_info[0])
if sys.version_info[0] < 3:
    print("executing python 2 imports")
    from BaseHTTPServer import HTTPServer
    from SimpleHTTPServer import SimpleHTTPRequestHandler
    from SocketServer import ThreadingMixIn
else:
    print("executing python 3 imports")
    from http.server import HTTPServer, SimpleHTTPRequestHandler
    from socketserver import ThreadingMixIn

import cgi
import os
import sys
import xml.etree.ElementTree as ET
import hashlib
import binascii


def convertToBits(msg):
    if sys.version[0] < 3:
        return msg
    else:
        return bytes(msg, "utf-8")


class ServerHandler(SimpleHTTPRequestHandler):

    def authenticate(self):
        method, authString = self.headers.getheader('Authorization').split(" ", 1)
        username, password = authString.decode("base64").split(":", 1)

        # parse accounts xml file
        global accounts_path
        filePath = accounts_path
        if '"' in filePath:
            filePath = filePath.replace('"', '')
        tree = ET.parse(filePath)
        root = tree.getroot()
        usernameEntry = None
        for account in root:
            usernameEntry = account.attrib
            if username == usernameEntry["username"]:
                # check the password
                salt = account[0].text
                storedpassword = account[1].text
                hashedPassword = binascii.hexlify(hashlib.pbkdf2_hmac('sha256', b'%s' % storedpassword,
                                                                    b'%s' % salt, 10000))
                if password == hashedPassword:
                    return True
                else:
                    return False
        return False

    def do_HEAD(self):
        self.send_response(200)
        self.send_header('Content-Type', 'text/html')
        self.end_headers()

    def do_AUTHHEAD(self):
        print("send header")
        self.send_response(401)
        self.send_header('WWW-Authenticate', 'Basic realm=\"Login to fileserver\"')
        self.send_header('Content-Type', 'text/html')
        self.end_headers()

    def do_GET(self):
        if self.headers.getheader('Authorization') == None:
            self.do_AUTHHEAD()
            self.wfile.write('no auth header received')
            pass
        else:
            # if os.environ.get("DBFILESERVER_ACCOUNTS_FILE") is None:
            #     self.send_response(500)
            #     self.end_headers()
            #     self.wfile.write("Cannot login, system error. No account information to check against.")
            if self.authenticate():
                SimpleHTTPRequestHandler.do_GET(self)
            else:
                self.do_AUTHHEAD()
                self.wfile.write(self.headers.getheader('Authorization'))
                self.wfile.write('not authenticated')
            pass

    def do_POST(self):
        # Parse the form data posted
        # if os.environ.get("DBFILESERVER_ACCOUNTS_FILE") is None:
        #         self.send_response(500)
        #         self.end_headers()
        #         self.wfile.write("Cannot login, system error. No account information to check against.")
        if not self.authenticate():
            self.do_AUTHHEAD()
            return
        form = cgi.FieldStorage(
            fp=self.rfile, 
            headers=self.headers,
            environ={'REQUEST_METHOD':'POST',
                     'CONTENT_TYPE':self.headers['Content-Type'],
                     })

        exception = None
        # write file to memory
        for field in form.keys():
            field_item = form[field]
            if field_item.filename:
                try:
                    filePath = os.path.join(os.getcwd(),
                                            self.path[1:],
                                            field_item.filename)

                    with open(filePath, "wb") as newFile:
                        newFile.write(field_item.file.read())
                except Exception as e:
                    print("Exception: %s" % e)
                    exception = e

        if exception is None:
            # Begin the response
            self.send_response(200)
            self.end_headers()
            self.wfile.write('Client: %s\n' % str(self.client_address))
            self.wfile.write('User-agent: %s\n' % str(self.headers['user-agent']))
            self.wfile.write('Path: %s\n' % self.path)
            self.wfile.write('Form data:\n')

            # Echo back information about what was posted in the form
            for field in form.keys():
                field_item = form[field]
                if field_item.filename:
                    # The field contains an uploaded file
                    file_data = field_item.file.read()
                    file_len = len(file_data)
                    del file_data
                    self.wfile.write('\tUploaded %s as "%s" (%d bytes)\n' % \
                            (field, field_item.filename, file_len))

                else:
                    # Regular form value
                    self.wfile.write('\t%s=%s\n' % (field, form[field].value))
        else:
            # send back error response
            self.send_response(500)
            self.end_headers()
            self.wfile.write("Error: %s" % exception)

        return

    def do_DELETE(self):
        # if os.environ.get("DBFILESERVER_ACCOUNTS_FILE") is None:
        #         self.send_response(500)
        #         self.end_headers()
        #         self.wfile.write("Cannot login, system error. No account information to check against.")
        if not self.authenticate():
            self.do_AUTHHEAD()
            return
        # Parse the form data posted
        filePath = os.path.join(os.getcwd(), self.path[1:])
        responseCode = 501
        try:
            if os.path.exists(filePath):
                os.remove(filePath)
                responseCode = 200
                msg = "File %s deleted (%s bytes)" % (filePath, numBytes)
            else:
                msg = "File %s does not exist" % filePath
        except Exception as e:
            msg = "Unknown error: %s" % e

        self.send_response(responseCode)
        self.end_headers()
        self.wfile.write(msg)
        return

    def do_LIST(self):
        # if os.environ.get("DBFILESERVER_ACCOUNTS_FILE") is None:
        #         self.send_response(500)
        #         self.end_headers()
        #         self.wfile.write("Cannot login, system error. No account information to check against.")
        if not self.authenticate():
            self.do_AUTHHEAD()
            return
        filePath = os.path.join(os.getcwd(), self.path[1:])
        responseCode = 200
        self.send_response(responseCode)
        self.end_headers()
        dirContents = "\n".join(os.listdir(filePath))
        self.wfile.write(bytes(dirContents))
        return

class MultithreadedServer(ThreadingMixIn,
                          HTTPServer):
    pass


def ParseCommandLine(args):
    commands = {}
    for arg in args:
        if arg.startswith("-"):
            vals = arg.split("=")
            if len(vals) == 2:
                commands[vals[0][1:]] = vals[1]
            elif len(vals) == 1:
                commands[vals[0][1:]] = True
    return commands


if __name__ == "__main__":
    commands = ParseCommandLine(sys.argv[1:])
    print("commands: %s" % commands)
    PORT = 8000
    IP = ""
    if "port" in commands:
        PORT = int(commands["port"])
    if "ip" in commands:
        IP = commands["ip"]
    if "wd" in commands:
        path = commands["wd"]
        if os.path.isabs(path):
            os.chdir(path)
        else:
            os.chdir(os.path.join(os.getcwd(), commands["wd"]))
    accounts_path = os.environ["DBFILESERVER_ACCOUNTS_PATH"]\
        if "config" not in commands or not\
           os.path.isabs(commands["config"]) else commands["config"]

    try:
        # this is NOT an instantiation, but I am storing the
        # type ServerHandler. MultithreadedServer, if given
        # the handler type, will instantiate it.
        Handler = ServerHandler

        server = MultithreadedServer(("", PORT), Handler)

        print("Serving at: http://%(interface)s:%(port)s" %
              dict(interface=IP or "localhost", port=PORT))
        server.serve_forever()
    except KeyboardInterrupt:
        print("shutting down file server")
        server.socket.close()
