# SYSTEM IMPORTS
import os
import requests
import sys

# PYTHON PROJECT IMPORTS

# takes components of a url path and converts them into
# a legit url
def urljoin(url, *urls):
    urlList = [url]
    urlList.extend([urlPart for urlPart in urls])
    unrefinedUrl = '/'.join(urlList).strip()
    unrefinedUrl = unrefinedUrl.replace("//", "/")
    return unrefinedUrl.replace("http:/", "http://")


class HTTPRequest(object):

    # baseUrl here is the url of the fileserver you are trying to connect to.
    def __init__(self, baseUrl):
        self.user = os.environ.get("DBFILESERVER_USERNAME")
        self.pswrd = os.environ.get("DBFILESERVER_PASSWORD")
        self.baseUrl = baseUrl

    # used to download a file from the fileserver
    # receivingFilePath is the file path you want to save the downloaded file as.
    # urlParams are the list of url parts that will point to the file you want to download
    #	 for example: if self.baseUrl = http://123.45:80/ and urlParams = ["foo", "myfile.txt"]
    # 				  then you would be downloading http://123.45:80/foo/myfile.txt
    # don't worry about the other arguments, you want to keep them set at default.
    def download(self, receivingFilePath, urlParams=[], fileChunkSize=1, readBytes=True):
        url = None
        if len(urlParams) == 0:
            url = self.baseUrl
        else:
            url = urljoin(self.baseUrl, *urlParams)
        response = requests.get(url, stream=True, auth=requests.auth.HTTPBasicAuth(self.user, self.pswrd))

        if response.status_code != 200:
            print("ERROR downoading file %s: %s" % (url, response))
            return

        numBytes = len(response.content)
        currentBytes = 0.0
        minPercentToPrint = 0
        print("Starting download (%s bytes):" % numBytes)
        with open(receivingFilePath, ("wb" if readBytes else "w")) as f:
            for chunk in response.iter_content(fileChunkSize):
                if currentBytes/numBytes >= minPercentToPrint:
                    if sys.version_info[0] < 3:
                        print("[%s%%]" % int(currentBytes/numBytes * 100)),
                    else:
                        print("[%s%%]" % int(currentBytes/numBytes * 100), end=" ")
                    minPercentToPrint += 0.1
                f.write(chunk)
                currentBytes += len(chunk)
        print("Download done")

    # upload a file to the fileserver
    # this works by providing a file path (to the file you want to upload)
    # and the url where you want to upload to. NOTE: you don't create the
    # name of the uploaded file, the file when uploaded will be the same
    # as the one you wanted to upload.
    # EXAMPLE:
    #		upload("C:\\Users\myFile.txt", urlParams=["myDirOnFileServer"])
    #		will upload "myFile.txt" to the fileserver at the url
    #		"http://123.45:80/myDirOnFileServer/myFile.txt"
    #		assuming self.baseUrl = "http://123.45:80/"
    def upload(self, filePath, fileName="", urlParams=[]):
        fullFilePath = None
        url = None
        if fileName == "":
            fullFilePath = filePath

            # parse out fileName from filePath
            fileName = filePath.split(os.sep)[-1]
        else:
            fullFilePath = os.path.join(filePath, fileName)

        if len(urlParams) == 0:
            url = self.baseUrl
        else:
            url = urljoin(self.baseUrl, *urlParams)

        # to post, do I have to add "/post" to the end of the url?
        response = requests.post(url, files={"upload_file": open(fullFilePath, 'rb')},
                                 auth=requests.auth.HTTPBasicAuth(self.user, self.pswrd))

        # handle response
        if response.status_code != 200:
            print("Error %s uploading %s to %s" % (response.status_code,
                                                   fullFilePath, url))

    # delete a file from the file server.
    # you specify the url you want to delete the same way as downloading a file
    def delete(self, urlParams=[]):
        fullUrlPath = self.baseUrl if len(urlParams) == 0 else urljoin(self.baseUrl, *urlParams)

        response = requests.delete(fullUrlPath, auth=requests.auth.HTTPBasicAuth(self.user, self.pswrd))
        if response.status_code != 200:
            print("Error %s deleting file at url %s" % (response.status_code,
                                                        fullUrlPath))

    # delete the contents of a directory.
    # you specify the directory url the same way as uploading a file.
    # NOTE this does not delete directories recursively...I'll work on that.
    def deleteAll(self, urlParams=[]):
        fullUrlPath = self.baseUrl if len(urlParams) == 0 else urljoin(self.baseUrl, *urlParams)

        for fileToDelete in self.listUrlContents(fullUrlPath).split("\n"):
            self.delete(urlParams + [fileToDelete])

    # get the contents of a url directory (a directory on the fileserver).
    # specify the directory just like uploading a file.
    def listUrlContents(self, urlParams=[]):
        fullUrlPath = self.baseUrl if len(urlParams) == 0 else urljoin(self.baseUrl, *urlParams)

        response = requests.request("LIST", fullUrlPath, auth=requests.auth.HTTPBasicAuth(self.user, self.pswrd))
        if response.status_code != 200:
            print("Error %s listing contents of %s" % (response.status_code,
                                                       fullUrlPath))
        if sys.version_info[0] < 3:
            return str(response.content)
        else:
            return response.content.decode('utf-8')
